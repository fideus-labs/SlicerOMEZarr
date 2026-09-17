"""OME-Zarr (OME-NGFF) support for 3D Slicer.

Registers:

* ``OMEZarrFileReader`` - scripted file reader (Add Data, drag-and-drop,
  ``slicer.util.loadNodeFromFile``), for local stores, ``.ozx`` files and URLs.
* ``OMEZarrFileWriter`` - scripted file writer for scalar and label map volumes.
* ``OMEZarrFileDialog`` - drop target for OME-Zarr directories.
* ``OMEZarrLogic`` - level selection, NGFF -> RAS geometry, chunk-limited reads,
  labels, time series, view-driven refinement.

NGFF parsing, multiscales, store access and RFC-4 orientation come from ngff-zarr.
"""

import json
import logging
import os
import re
import urllib.request

import numpy as np
import qt
import vtk

import slicer
import slicer.packaging
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)

SPATIAL_DIMS = ("x", "y", "z")

# Slicer world coordinates are millimetres. NGFF axis units are UDUNITS-2 names.
LENGTH_UNIT_TO_MM = {
    None: 1.0,
    "": 1.0,
    "millimeter": 1.0,
    "micrometer": 1e-3,
    "micron": 1e-3,
    "nanometer": 1e-6,
    "picometer": 1e-9,
    "angstrom": 1e-7,
    "centimeter": 10.0,
    "decimeter": 100.0,
    "meter": 1000.0,
    "kilometer": 1e6,
    "inch": 25.4,
    "foot": 304.8,
    "yard": 914.4,
}

# Symbols shown in the panel, and used when the display follows the store's unit.
UNIT_SYMBOLS = {"micrometer": "µm", "micron": "µm", "nanometer": "nm", "millimeter": "mm", "centimeter": "cm"}

DEFAULT_BUDGET_FRACTION = 0.25  # of available RAM when the budget setting is "auto"
FALLBACK_MAX_BYTES = 1 << 30

# Slicer core lookup tables used to colour separate microscopy channels.
CHANNEL_COLOR_NODE_IDS = {
    "red": "vtkMRMLColorTableNodeRed",
    "green": "vtkMRMLColorTableNodeGreen",
    "blue": "vtkMRMLColorTableNodeBlue",
    "yellow": "vtkMRMLColorTableNodeYellow",
    "cyan": "vtkMRMLColorTableNodeCyan",
    "magenta": "vtkMRMLColorTableNodeMagenta",
    "grey": "vtkMRMLColorTableNodeGrey",
}

MAX_COLOR_TABLE_ENTRIES = 65536


class Settings:
    """QSettings keys. Values are read on each use so the module panel takes effect immediately."""

    MAX_BYTES = "OMEZarr/MaxBytes"  # int bytes, 0 = automatic
    ORIENTATION = "OMEZarr/AssumedOrientation"  # LPS or RAS, used when a store has no RFC-4 orientation
    LOAD_LABELS = "OMEZarr/LoadLabels"
    TIME_MODE = "OMEZarr/TimeMode"  # sequence or index
    DISPLAY_UNITS = "OMEZarr/DisplayUnits"  # switch Slicer's length display unit to the store's
    AUTO_REFINE = "OMEZarr/AutoRefine"  # refine the slice views automatically after a store is loaded
    LABELS_AS_SEGMENTATION = "OMEZarr/LabelsAsSegmentation"  # load labels as Segmentation nodes
    STORAGE_OPTIONS = "OMEZarr/StorageOptions"  # JSON passed to ngff-zarr for remote stores
    DETECT_LABEL_MAPS = "OMEZarr/DetectLabelMaps"  # integer stores with few values load as label maps

    @staticmethod
    def get(key, default):
        value = qt.QSettings().value(key)
        if value is None or value == "":
            return default
        if isinstance(default, bool):
            return str(value).lower() in ("true", "1", "yes")
        if isinstance(default, int):
            try:
                return int(value)
            except (TypeError, ValueError):
                return default
        return str(value)

    @staticmethod
    def set(key, value):
        qt.QSettings().setValue(key, value)


#
# Module
#


class OMEZarr(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("OME-Zarr")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "Informatics")]
        self.parent.dependencies = ["Sequences"]
        self.parent.contributors = ["Valentin Boussot (Fideus Labs)", "Matt McCormick (Fideus Labs)"]
        self.parent.helpText = _(
            "Read and write OME-Zarr (OME-NGFF) images. Drag an .ome.zarr directory into Slicer, or use "
            "File > Add Data. The finest resolution level fitting the memory budget is loaded; "
            "labels and time series are loaded alongside. Use 'Refine current view' to reload what a "
            "slice view shows at a finer level."
        )
        self.parent.acknowledgementText = _("Built on ngff-zarr (https://github.com/fideus-labs/ngff-zarr).")


#
# Store discovery helpers
#


def isRemoteUrl(path):
    return bool(re.match(r"^(https?|s3|gs|gcs|az|abfs)://", str(path)))


def normalizeStorePath(path):
    """Canonical form of a store path for node attributes and comparisons (URLs unchanged)."""
    path = str(path)
    if isRemoteUrl(path):
        return path.rstrip("/")
    return os.path.normpath(os.path.abspath(path))


def samePath(a, b):
    return a is not None and b is not None and normalizeStorePath(a) == normalizeStorePath(b)


def joinStorePath(root, *parts):
    root = str(root).rstrip("/")
    return "/".join([root, *parts]) if isRemoteUrl(root) else os.path.join(root, *parts)


def readStoreAttributes(root):
    """Merged root attributes of a Zarr group (v2 ``.zattrs`` or v3 ``zarr.json``), OME keys unwrapped.

    Returns None when the group has no readable metadata. Local paths and http(s) URLs only.
    """
    root = str(root).rstrip("/")
    for name in ("zarr.json", ".zattrs"):
        target = joinStorePath(root, name)
        try:
            if isRemoteUrl(target):
                if not target.startswith("http"):
                    return None
                with urllib.request.urlopen(target, timeout=20) as response:  # noqa: S310 - user-provided URL
                    attrs = json.loads(response.read().decode("utf-8"))
            else:
                if not os.path.isfile(target):
                    continue
                with open(target, encoding="utf-8") as fp:
                    attrs = json.load(fp)
        except (OSError, ValueError):
            continue
        if not isinstance(attrs, dict):
            continue
        if name == "zarr.json":
            attrs = attrs.get("attributes", {})
        merged = dict(attrs)
        if isinstance(attrs.get("ome"), dict):
            merged.update(attrs["ome"])
        return merged
    return None


def omeZarrRootFromPath(path):
    """Return the OME-Zarr multiscales root for a path, or None.

    Accepts the store directory itself, one of its metadata files
    (``zarr.json``/``.zattrs``, which is what ``Add Data`` lists when a directory
    is added), a ``.ozx``/``.zip`` file, or a remote URL.
    """
    if not path:
        return None
    path = str(path)
    if isRemoteUrl(path):
        return path.rstrip("/") if re.search(r"\.zarr/?$", path) else None
    path = os.path.abspath(path)
    if os.path.isfile(path):
        if os.path.basename(path) in ("zarr.json", ".zattrs", ".zgroup"):
            path = os.path.dirname(path)
        elif path.lower().endswith((".ozx", ".zarr.zip")):
            return path
        else:
            return None
    if not os.path.isdir(path):
        return None
    attrs = readStoreAttributes(path)
    return path if attrs and ("multiscales" in attrs or "bioformats2raw.layout" in attrs) else None


def isBioformats2rawRoot(root):
    attrs = readStoreAttributes(root)
    return bool(attrs and "bioformats2raw.layout" in attrs and "multiscales" not in attrs)


def bioformats2rawSeries(root, limit=1000):
    """Paths of the image series of a bioformats2raw container: ``<root>/0``, ``<root>/1``, ..."""
    series = []
    for index in range(limit):
        candidate = joinStorePath(root, str(index))
        attrs = readStoreAttributes(candidate)
        if not attrs or "multiscales" not in attrs:
            break
        series.append(candidate)
    return series


def isLabelStore(root):
    attrs = readStoreAttributes(root)
    return bool(attrs and "image-label" in attrs)


def labelGroupNames(root):
    """Names listed in the ``labels`` group of an image store (empty when absent)."""
    attrs = readStoreAttributes(joinStorePath(root, "labels"))
    names = attrs.get("labels") if attrs else None
    return [str(n) for n in names] if isinstance(names, list) else []


def scalarVolumes(nodes):
    return [n for n in nodes if n.IsA("vtkMRMLScalarVolumeNode") and not n.IsA("vtkMRMLLabelMapVolumeNode")]


def labelMaps(nodes):
    return [n for n in nodes if n.IsA("vtkMRMLLabelMapVolumeNode")]


def runResponsive(work, onTick=None):
    """Run ``work()`` in a worker thread while this (GUI) thread keeps processing Qt events.

    ``onTick`` is called from the GUI thread about a hundred times a second. An exception
    raised by ``work`` is re-raised here; its return value is returned.
    """
    import threading
    import time

    outcome = {}

    def target():
        try:
            outcome["result"] = work()
        except BaseException as e:  # noqa: BLE001 - re-raised in the calling thread
            outcome["error"] = e

    thread = threading.Thread(target=target, name="OMEZarr", daemon=True)
    thread.start()
    while thread.is_alive():
        slicer.app.processEvents()
        if onTick:
            onTick()
        time.sleep(0.01)
    thread.join()
    if "error" in outcome:
        raise outcome["error"]
    return outcome.get("result")


#
# Progress reporting
#


class Progress:
    """Progress dialog with a cancel button when the GUI is up; silent otherwise.

    Used as ``progress(done, total, text)``; returns False once the user cancelled.
    """

    def __init__(self, label):
        self.label = label
        self.dialog = None

    def __enter__(self):
        if slicer.util.mainWindow() and not slicer.app.testingEnabled():
            self.dialog = slicer.util.createProgressDialog(labelText=self.label, windowTitle=_("OME-Zarr"), maximum=100)
        return self

    def __call__(self, done, total, text=None):
        if self.dialog is None:
            return True
        self.dialog.value = int(100 * done / max(1, total))
        if text:
            self.dialog.labelText = text
        return not self.dialog.wasCanceled

    def __exit__(self, *args):
        if self.dialog is not None:
            self.dialog.close()
        return False


#
# Logic
#


class OMEZarrLogic(ScriptedLoadableModuleLogic):
    """NGFF <-> MRML mapping. Stateless apart from a small multiscales cache."""

    _multiscalesCache = {}

    @staticmethod
    def ensureNgffZarr():
        # From 0.46.1 on, ngff-zarr fills NgffImage.axes_orientations from the RFC-4 metadata on read.
        requirement = "ngff-zarr[remote]>=0.46.1"
        if not slicer.packaging.pip_check(requirement):
            interactive = slicer.util.mainWindow() and not slicer.app.testingEnabled()
            if interactive and not slicer.util.confirmOkCancelDisplay(
                _("ngff-zarr 0.46.1 or newer is required to read OME-Zarr images. Install it now?")
            ):
                raise RuntimeError("ngff-zarr 0.46.1 or newer is not installed")
            with slicer.util.tryWithErrorDisplay(_("Failed to install ngff-zarr"), waitCursor=True):
                slicer.util.pip_install(requirement)
        import ngff_zarr

        return ngff_zarr

    @staticmethod
    def storageOptions(path):
        """Options for a remote store: the configured JSON, and anonymous S3 access when no
        AWS credentials are configured (otherwise obstore spends half a minute asking the EC2
        metadata service before failing on a public bucket)."""
        if not isRemoteUrl(path):
            return None
        options = {}
        configured = Settings.get(Settings.STORAGE_OPTIONS, "")
        if configured:
            try:
                options.update(json.loads(configured))
            except ValueError:
                logging.warning("Ignoring invalid JSON in the OME-Zarr storage options setting")
        if str(path).startswith("s3://") and "anon" not in options and "skip_signature" not in options:
            hasCredentials = any(
                os.environ.get(name) for name in ("AWS_ACCESS_KEY_ID", "AWS_PROFILE", "AWS_SESSION_TOKEN")
            ) or os.path.isfile(os.path.expanduser("~/.aws/credentials"))
            if not hasCredentials:
                options["anon"] = True
        return options or None

    @classmethod
    def openMultiscales(cls, path, useCache=True):
        """Open a store with ngff-zarr; arrays stay lazy (dask), nothing is read yet."""
        ngff_zarr = cls.ensureNgffZarr()
        key = str(path)
        if useCache and key in cls._multiscalesCache:
            return cls._multiscalesCache[key]
        source = key
        if isBioformats2rawRoot(key):
            series = bioformats2rawSeries(key)
            if not series:
                raise ValueError(f"No image series found in the bioformats2raw container {key}")
            source = series[0]
        options = cls.storageOptions(source)
        multiscales = (
            ngff_zarr.from_ome_zarr(source, storage_options=options) if options else ngff_zarr.from_ome_zarr(source)
        )
        if useCache:
            cls._multiscalesCache[key] = multiscales
        return multiscales

    @classmethod
    def clearCache(cls):
        cls._multiscalesCache.clear()

    # ---- metadata ----

    @staticmethod
    def spatialShape(image):
        dims = list(image.dims)
        return {d: image.data.shape[dims.index(d)] for d in SPATIAL_DIMS if d in dims}

    @staticmethod
    def axisLength(image, dim):
        dims = list(image.dims)
        return int(image.data.shape[dims.index(dim)]) if dim in dims else 1

    @classmethod
    def volumeBytes(cls, image, region=None):
        """Bytes of one spatial volume (single channel, single time point), or of ``region`` within it."""
        size = 1
        for d, n in cls.spatialShape(image).items():
            if region and d in region:
                n = region[d][1] - region[d][0]
            size *= int(n)
        return size * int(np.dtype(image.data.dtype).itemsize)

    @classmethod
    def levelInfo(cls, multiscales):
        info = []
        base = multiscales.images[0]
        for index, image in enumerate(multiscales.images):
            downsample = {d: image.scale.get(d, 1.0) / base.scale.get(d, 1.0) for d in SPATIAL_DIMS if d in image.dims}
            info.append(
                {
                    "level": index,
                    "dims": tuple(image.dims),
                    "shape": tuple(int(n) for n in image.data.shape),
                    "chunks": tuple(int(n) for n in image.data.chunksize),
                    "dtype": str(image.data.dtype),
                    "bytes": cls.volumeBytes(image),
                    "scale": dict(image.scale),
                    "units": dict(image.axes_units or {}),
                    "downsample": downsample,
                }
            )
        return info

    @classmethod
    def selectLevel(cls, multiscales, maxBytes, copies=1):
        """Finest level whose volume (times ``copies``: channels, time points) fits the budget."""
        for index, image in enumerate(multiscales.images):
            if cls.volumeBytes(image) * copies <= maxBytes:
                return index
        return len(multiscales.images) - 1

    @staticmethod
    def maxBytesFromSettings():
        """Configured budget, or a fraction of the available RAM when set to automatic."""
        value = Settings.get(Settings.MAX_BYTES, 0)
        if value > 0:
            return value
        try:
            import psutil

            return int(psutil.virtual_memory().available * DEFAULT_BUDGET_FRACTION)
        except Exception:  # noqa: BLE001 - psutil is optional
            return FALLBACK_MAX_BYTES

    @staticmethod
    def timeUnit(image):
        return str((image.axes_units or {}).get("t") or "")

    @staticmethod
    def lengthUnit(image):
        units = image.axes_units or {}
        return next((str(units[d]) for d in SPATIAL_DIMS if units.get(d)), "")

    # ---- geometry ----

    @staticmethod
    def unitScaleToMm(image, userMessages=None):
        factors = {}
        for d in SPATIAL_DIMS:
            unit = (image.axes_units or {}).get(d)
            unit = str(unit) if unit is not None else None
            if unit in LENGTH_UNIT_TO_MM:
                factors[d] = LENGTH_UNIT_TO_MM[unit]
            else:
                factors[d] = 1.0
                message = f"Unknown length unit '{unit}' for axis '{d}', values are used as millimetres."
                logging.warning(message)
                if userMessages:
                    userMessages.AddMessage(vtk.vtkCommand.WarningEvent, message)
        return factors

    @classmethod
    def ijkToRasMatrix(cls, image, userMessages=None):
        """4x4 IJK->RAS (mm) and the orientation source ("rfc4", "assumed-LPS", "assumed-RAS").

        NGFF x/y/z are ITK/LPS physical axes (the ngff-zarr convention) and ``translation``
        is the origin. RFC-4 orientation of every spatial axis gives the direction cosines.
        Without it, the axes are assumed LPS or RAS per the module setting.
        """
        from ngff_zarr.rfc4 import anatomical_orientation_to_itk_direction

        dims = list(image.dims)
        factors = cls.unitScaleToMm(image, userMessages)
        spatialDims = [d for d in SPATIAL_DIMS if d in dims]
        orientations = image.axes_orientations or {}
        columns = {}
        if spatialDims and all(d in orientations for d in spatialDims):
            for d in spatialDims:
                column = anatomical_orientation_to_itk_direction(orientations[d].value)
                if column is None:
                    columns = {}
                    break
                columns[d] = column
        spacing = [image.scale.get(d, 1.0) * factors[d] if d in dims else 1.0 for d in SPATIAL_DIMS]
        origin = [image.translation.get(d, 0.0) * factors[d] if d in dims else 0.0 for d in SPATIAL_DIMS]
        direction = np.eye(3)
        if columns:
            source = "rfc4"
            for columnIndex, d in enumerate(SPATIAL_DIMS):
                if d in columns:
                    direction[:, columnIndex] = columns[d]
            toRas = np.diag([-1.0, -1.0, 1.0, 1.0])
        else:
            assumed = Settings.get(Settings.ORIENTATION, "LPS").upper()
            source = f"assumed-{assumed}"
            toRas = np.eye(4) if assumed == "RAS" else np.diag([-1.0, -1.0, 1.0, 1.0])
            if userMessages:
                userMessages.AddMessage(
                    vtk.vtkCommand.MessageEvent,
                    f"No anatomical orientation (RFC-4) in the store; x/y/z axes assumed {assumed} "
                    "(OME-Zarr module settings).",
                )
        ijkToPhysical = np.eye(4)
        ijkToPhysical[:3, :3] = direction @ np.diag(spacing)
        ijkToPhysical[:3, 3] = origin
        return toRas @ ijkToPhysical, source

    @classmethod
    def regionIjkToRas(cls, image, region=None, userMessages=None):
        """IJK->RAS of ``region`` within the image (of the whole image without a region), and its source."""
        ijkToRas, source = cls.ijkToRasMatrix(image, userMessages)
        if region:
            start = np.array([region.get(d, (0, None))[0] for d in SPATIAL_DIMS], dtype=float)
            ijkToRas[:3, 3] = (ijkToRas @ np.append(start, 1.0))[:3]
        return ijkToRas, source

    # ---- array access ----

    @staticmethod
    def computeArray(darray, output, progress=None, label=""):
        """Compute a dask array into ``output`` without freezing the application.

        The array is read slab by slab along its first axis (groups of chunks of about
        64 MiB, at most 32), each slab a parallel dask compute, in a worker thread.
        ``output`` is preallocated, typically a view on a vtkImageData buffer, so a volume is
        never held twice in memory. Cancelling stops after the slab being read and raises
        InterruptedError.
        """
        chunks = list(darray.chunks[0])
        steps = max(1, min(32, int(np.ceil(darray.nbytes / (64 << 20)))))
        groupSize = max(1, int(np.ceil(len(chunks) / steps)))
        total = int(np.ceil(len(chunks) / groupSize))
        state = {"done": 0, "cancel": False}

        def work():
            start = 0
            for step in range(total):
                if state["cancel"]:
                    return
                size = sum(chunks[step * groupSize : (step + 1) * groupSize])
                output[start : start + size] = darray[start : start + size].compute()
                start += size
                state["done"] = step + 1

        def onTick():
            if progress and not progress(state["done"], total, label):
                state["cancel"] = True

        runResponsive(work, onTick)
        if state["cancel"]:
            raise InterruptedError("Loading cancelled")

    @classmethod
    def spatialDaskArray(cls, image, timeIndex=0, channelIndex=0, region=None, userMessages=None):
        """The lazy (z, y, x) sub-array of one time point and channel, and whether a z axis must be added.

        ``region`` maps spatial dim -> (start, stop) in this level's index space; only the
        chunks it intersects are read when the array is computed.
        """
        dims = list(image.dims)
        index = []
        for d in dims:
            if d == "t":
                index.append(int(timeIndex))
            elif d == "c":
                index.append(int(channelIndex))
            elif d in SPATIAL_DIMS:
                index.append(slice(*region[d]) if region and d in region else slice(None))
            else:
                message = f"Axis '{d}' is not supported, using its first element."
                logging.warning(message)
                if userMessages:
                    userMessages.AddMessage(vtk.vtkCommand.WarningEvent, message)
                index.append(0)
        sub = image.data[tuple(index)]
        remaining = [d for d in dims if d in SPATIAL_DIMS]
        if "x" not in remaining or "y" not in remaining:
            raise ValueError("OME-Zarr image must have x and y axes")
        order = [remaining.index(d) for d in ("z", "y", "x") if d in remaining]
        sub = sub.transpose(order)
        return sub, "z" not in remaining

    @staticmethod
    def vtkCompatibleDtype(dtype):
        dtype = np.dtype(dtype)
        if dtype == np.bool_:
            return np.dtype(np.uint8)
        if dtype == np.float16:
            return np.dtype(np.float32)
        return dtype

    @classmethod
    def fillVolumeNode(
        cls, node, image, timeIndex, channelIndex, region, ijkToRas, userMessages, progress, label, dtype=None
    ):
        """Read a (z, y, x) volume straight into the node's image buffer.

        The vtkImageData is allocated first and the dask array is computed into a numpy
        view of its scalars, so peak memory is one copy of the volume, not two. ``dtype``
        overrides the stored type (label maps must be integer).
        """
        from vtk.util import numpy_support

        sub, addZ = cls.spatialDaskArray(image, timeIndex, channelIndex, region, userMessages)
        shape = (1, *sub.shape) if addZ else tuple(sub.shape)
        dtype = np.dtype(dtype) if dtype else cls.vtkCompatibleDtype(sub.dtype)
        imageData = vtk.vtkImageData()
        imageData.SetDimensions(int(shape[2]), int(shape[1]), int(shape[0]))
        imageData.AllocateScalars(numpy_support.get_vtk_array_type(dtype), 1)
        view = numpy_support.vtk_to_numpy(imageData.GetPointData().GetScalars()).reshape(shape)
        cls.computeArray(sub.astype(dtype), view[0] if addZ else view, progress, label)
        node.SetAndObserveImageData(imageData)
        node.SetIJKToRASMatrix(slicer.util.vtkMatrixFromArray(ijkToRas))
        return node

    # ---- channels and display ----

    @staticmethod
    def channelDescriptions(multiscales, image):
        """List of {label, color, window} per channel, from OMERO metadata when present."""
        dims = list(image.dims)
        count = int(image.data.shape[dims.index("c")]) if "c" in dims else 1
        descriptions = [{"label": None, "color": None, "window": None} for _ in range(count)]
        names = list(image.channel_names or [])
        colors = list(image.channel_colors or [])
        for i in range(count):
            if i < len(names) and names[i]:
                descriptions[i]["label"] = str(names[i])
            if i < len(colors) and colors[i]:
                descriptions[i]["color"] = str(colors[i])
        omero = getattr(multiscales.metadata, "omero", None)
        if omero is not None and getattr(omero, "channels", None):
            for i, channel in enumerate(omero.channels[:count]):
                if getattr(channel, "label", None):
                    descriptions[i]["label"] = channel.label
                if getattr(channel, "color", None):
                    descriptions[i]["color"] = channel.color
                window = getattr(channel, "window", None)
                if window is not None:
                    descriptions[i]["window"] = (float(window.start), float(window.end))
        return descriptions

    @staticmethod
    def colorNodeIdForHexColor(hexColor):
        try:
            r, g, b = (int(hexColor.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
        except (ValueError, AttributeError, TypeError):
            return None
        high = {name for name, value in (("r", r), ("g", g), ("b", b)) if value >= 128}
        return CHANNEL_COLOR_NODE_IDS.get(
            {
                frozenset("r"): "red",
                frozenset("g"): "green",
                frozenset("b"): "blue",
                frozenset("rg"): "yellow",
                frozenset("gb"): "cyan",
                frozenset("rb"): "magenta",
                frozenset("rgb"): "grey",
            }.get(frozenset(high))
        )

    @staticmethod
    def setupDisplay(node, description, multiChannel):
        node.CreateDefaultDisplayNodes()
        display = node.GetDisplayNode()
        if display is None:
            return
        if description.get("window"):
            start, end = description["window"]
            if end > start:
                display.SetAutoWindowLevel(False)
                display.SetWindowLevelMinMax(start, end)
        colorNodeId = OMEZarrLogic.colorNodeIdForHexColor(description.get("color")) if multiChannel else None
        if colorNodeId and slicer.mrmlScene.GetNodeByID(colorNodeId):
            display.SetAndObserveColorNodeID(colorNodeId)

    @staticmethod
    def defaultNodeName(path):
        base = os.path.basename(str(path).rstrip("/"))
        for suffix in (".ome.zarr", ".zarr.zip", ".ozx", ".zarr"):
            if base.lower().endswith(suffix):
                return base[: -len(suffix)]
        return base or "OMEZarr"

    @staticmethod
    def setNodeAttributes(node, path, level, dims, timeIndex, channelIndex, lengthUnit, orientationSource, region):
        node.SetAttribute("OMEZarr.Path", normalizeStorePath(path))
        node.SetAttribute("OMEZarr.Level", str(level))
        node.SetAttribute("OMEZarr.Dims", ",".join(dims))
        node.SetAttribute("OMEZarr.TimeIndex", str(timeIndex))
        node.SetAttribute("OMEZarr.Channel", str(channelIndex))
        node.SetAttribute("OMEZarr.LengthUnit", lengthUnit)
        node.SetAttribute("OMEZarr.OrientationSource", orientationSource)
        if region:
            node.SetAttribute("OMEZarr.Region", ";".join(f"{d}:{s}-{e}" for d, (s, e) in region.items()))

    # ---- loading ----

    @classmethod
    def loadImage(
        cls,
        path,
        level=None,
        timeIndex=None,
        channels=None,
        region=None,
        name=None,
        maxBytes=None,
        userMessages=None,
        multiscales=None,
        labels=None,
        timeMode=None,
        progress=None,
        asLabelMap=None,
    ):
        """Load a store into volume nodes. Returns the created nodes (volumes, proxies, label maps).

        ``region``: spatial dim -> (start, stop) at the chosen level (partial read).
        ``timeIndex``: a single time point; otherwise the setting decides between a
        Sequence of all time points and index 0. ``labels``: also load the ``labels`` groups.
        """
        if multiscales is None and isBioformats2rawRoot(path):
            series = bioformats2rawSeries(path)
            if not series:
                raise ValueError(f"No image series found in the bioformats2raw container {path}")
            baseName = name or cls.defaultNodeName(path)
            nodes = []
            for index, seriesPath in enumerate(series):
                seriesName = baseName if len(series) == 1 else f"{baseName}_{index}"
                nodes += cls.loadImage(
                    seriesPath,
                    level,
                    timeIndex,
                    channels,
                    region,
                    seriesName,
                    maxBytes,
                    userMessages,
                    None,
                    labels,
                    timeMode,
                    progress,
                    asLabelMap,
                )
            return nodes
        multiscales = multiscales or cls.openMultiscales(path)
        if isLabelStore(path) or (asLabelMap is None and cls.looksLikeLabelMap(multiscales)) or asLabelMap:
            return cls.loadLabelStore(
                path, level, region, name, userMessages, multiscales, progress=progress, maxBytes=maxBytes
            )
        baseImage = multiscales.images[0]
        descriptions = cls.channelDescriptions(multiscales, baseImage)
        channels = list(range(len(descriptions))) if channels is None else list(channels)
        timePoints = cls.axisLength(baseImage, "t")
        if timeIndex is None:
            timeMode = timeMode or Settings.get(Settings.TIME_MODE, "sequence")
            asSequence = timePoints > 1 and timeMode == "sequence"
            timeIndex = 0
        else:
            asSequence = False
            timeIndex = int(timeIndex)
        copies = len(channels) * (timePoints if asSequence else 1)
        if level is None:
            level = cls.selectLevel(multiscales, maxBytes or cls.maxBytesFromSettings(), copies)
        level = int(level)
        if level < 0 or level >= len(multiscales.images):
            raise ValueError(f"Level {level} out of range (0..{len(multiscales.images) - 1})")
        image = multiscales.images[level]
        dims = list(image.dims)

        if (
            region is None
            and userMessages
            and cls.volumeBytes(image) * copies > (maxBytes or cls.maxBytesFromSettings())
        ):
            userMessages.AddMessage(
                vtk.vtkCommand.WarningEvent,
                f"No resolution level fits the memory budget; loading level {level} "
                f"({cls.volumeBytes(image) * copies / 2**30:.2f} GiB). Use a region of interest for large stores.",
            )
        if level > 0 and region is None and userMessages:
            full = cls.volumeBytes(multiscales.images[0]) * copies
            factor = cls.levelInfo(multiscales)[level]["downsample"]
            userMessages.AddMessage(
                vtk.vtkCommand.WarningEvent,
                f"Loaded multiscale level {level} (downsampled by "
                f"{', '.join(f'{d}: {v:g}' for d, v in factor.items())}) because level 0 "
                f"would need {full / 2**30:.2f} GiB. Use 'Refine current view' or a region of "
                "interest in the OME-Zarr module for full resolution, or raise the memory budget.",
            )
        if timePoints > 1 and not asSequence and userMessages:
            userMessages.AddMessage(
                vtk.vtkCommand.MessageEvent, f"Time axis has {timePoints} points; loaded index {timeIndex}."
            )

        ijkToRas, orientationSource = cls.regionIjkToRas(image, region, userMessages)

        baseName = name or cls.defaultNodeName(path)
        lengthUnit = cls.lengthUnit(image)
        nodes = []
        sequences = []
        for channelIndex in channels:
            label = descriptions[channelIndex]["label"] or (f"c{channelIndex}" if len(descriptions) > 1 else None)
            nodeName = f"{baseName}_{label}" if label else baseName
            if region:
                nodeName += "_ROI"
            nodeName = slicer.mrmlScene.GenerateUniqueName(nodeName)
            if asSequence:
                sequence = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceNode", nodeName + "_sequence")
                sequence.SetIndexName("time")
                sequence.SetIndexUnit(cls.timeUnit(image))
                timeScale = float(image.scale.get("t", 1.0))
                for t in range(timePoints):
                    dataNode = slicer.vtkMRMLScalarVolumeNode()
                    dataNode.SetName(nodeName)
                    cls.fillVolumeNode(
                        dataNode,
                        image,
                        t,
                        channelIndex,
                        region,
                        ijkToRas,
                        userMessages,
                        progress,
                        f"{nodeName}, t={t + 1}/{timePoints}",
                    )
                    cls.setNodeAttributes(
                        dataNode, path, level, dims, t, channelIndex, lengthUnit, orientationSource, region
                    )
                    sequence.SetDataNodeAtValue(dataNode, f"{t * timeScale:g}")
                cls.setNodeAttributes(
                    sequence, path, level, dims, -1, channelIndex, lengthUnit, orientationSource, region
                )
                sequences.append((sequence, descriptions[channelIndex]))
            else:
                node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", nodeName)
                cls.fillVolumeNode(
                    node, image, timeIndex, channelIndex, region, ijkToRas, userMessages, progress, nodeName
                )
                cls.setNodeAttributes(
                    node, path, level, dims, timeIndex, channelIndex, lengthUnit, orientationSource, region
                )
                cls.setupDisplay(node, descriptions[channelIndex], multiChannel=len(descriptions) > 1)
                nodes.append(node)

        if sequences:
            browser = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceBrowserNode", baseName + "_browser")
            for index, (sequence, _description) in enumerate(sequences):
                if index == 0:
                    browser.SetAndObserveMasterSequenceNodeID(sequence.GetID())
                else:
                    browser.AddSynchronizedSequenceNodeID(sequence.GetID())
            slicer.modules.sequences.logic().UpdateProxyNodesFromSequences(browser)
            for sequence, description in sequences:
                proxy = browser.GetProxyNode(sequence)
                if proxy is not None:
                    cls.setupDisplay(proxy, description, multiChannel=len(descriptions) > 1)
                    nodes.append(proxy)

        if labels is None:
            labels = Settings.get(Settings.LOAD_LABELS, True)
        if labels:
            regionBounds = cls.rasBoundsOfRegion(image, region, ijkToRas) if region else None
            nodes += cls.loadLabelGroups(path, image, level, baseName, regionBounds, userMessages, progress)

        cls.applyDisplayUnits(lengthUnit)
        return nodes

    @staticmethod
    def rasBoundsOfRegion(image, region, ijkToRas):
        """RAS bounds of a region whose ``ijkToRas`` already has the region origin."""
        size = [region[d][1] - region[d][0] if d in region else OMEZarrLogic.axisLength(image, d) for d in SPATIAL_DIMS]
        corners = np.array(
            [[i, j, k, 1.0] for i in (0, size[0] - 1) for j in (0, size[1] - 1) for k in (0, size[2] - 1)]
        )
        ras = (ijkToRas @ corners.T).T[:, :3]
        return [ras[:, 0].min(), ras[:, 0].max(), ras[:, 1].min(), ras[:, 1].max(), ras[:, 2].min(), ras[:, 2].max()]

    # ---- labels ----

    @classmethod
    def matchingLevel(cls, multiscales, referenceImage):
        """Level of ``multiscales`` whose spacing is closest to ``referenceImage``'s."""
        reference = np.array([referenceImage.scale.get(d, 1.0) for d in SPATIAL_DIMS if d in referenceImage.dims])
        best, bestDistance = 0, None
        for index, image in enumerate(multiscales.images):
            scale = np.array([image.scale.get(d, 1.0) for d in SPATIAL_DIMS if d in referenceImage.dims])
            distance = (
                float(np.abs(np.log(scale / reference)).sum()) if scale.shape == reference.shape else float("inf")
            )
            if bestDistance is None or distance < bestDistance:
                best, bestDistance = index, distance
        return best

    @classmethod
    def loadLabelGroups(cls, root, image, level, baseName, rasBounds=None, userMessages=None, progress=None):
        nodes = []
        for labelName in labelGroupNames(root):
            labelPath = joinStorePath(root, "labels", labelName)
            try:
                labelMultiscales = cls.openMultiscales(labelPath)
                labelLevel = cls.matchingLevel(labelMultiscales, image)
                region = None
                if rasBounds is not None:
                    region = cls.regionFromRasBounds(labelMultiscales.images[labelLevel], rasBounds)
                nodes += cls.loadLabelStore(
                    labelPath, labelLevel, region, f"{baseName}_{labelName}", userMessages, labelMultiscales, progress
                )
            except Exception as e:  # noqa: BLE001 - a broken label group must not block the image
                message = f"Could not load labels '{labelName}': {e}"
                logging.exception(message)
                if userMessages:
                    userMessages.AddMessage(vtk.vtkCommand.WarningEvent, message)
        return nodes

    @classmethod
    def looksLikeLabelMap(cls, multiscales, maxValues=64):
        """Integer store without OMERO metadata whose coarsest level has few distinct values.

        Segmentation masks are often written as plain OME-Zarr images without
        ``image-label`` metadata; this heuristic loads them as label maps when the
        setting is on. The coarsest level is small, so the check is cheap.
        """
        if not Settings.get(Settings.DETECT_LABEL_MAPS, True):
            return False
        image = multiscales.images[-1]
        dtype = np.dtype(image.data.dtype)
        if not np.issubdtype(dtype, np.integer) or dtype.itemsize > 2 or "c" in image.dims or "t" in image.dims:
            return False
        if getattr(multiscales.metadata, "omero", None) is not None:
            return False
        if image.data.nbytes > (256 << 20):
            return False
        try:
            values = np.unique(np.asarray(image.data.compute()))
        except Exception:  # noqa: BLE001 - a failed probe just means "not a label map"
            return False
        return 1 < len(values) <= maxValues and values.min() >= 0

    @classmethod
    def loadLabelStore(
        cls,
        path,
        level=None,
        region=None,
        name=None,
        userMessages=None,
        multiscales=None,
        progress=None,
        maxBytes=None,
    ):
        """Load an ``image-label`` multiscales (or a store that looks like one) into a label map."""
        multiscales = multiscales or cls.openMultiscales(path)
        if level is None:
            level = cls.selectLevel(multiscales, maxBytes or cls.maxBytesFromSettings())
        image = multiscales.images[int(level)]
        dims = list(image.dims)
        ijkToRas, orientationSource = cls.regionIjkToRas(image, region, userMessages)
        nodeName = slicer.mrmlScene.GenerateUniqueName((name or cls.defaultNodeName(path)) + ("_ROI" if region else ""))
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", nodeName)
        stored = np.dtype(image.data.dtype)
        dtype = None if np.issubdtype(stored, np.integer) or stored == np.bool_ else np.int32
        cls.fillVolumeNode(node, image, 0, 0, region, ijkToRas, userMessages, progress, nodeName, dtype)
        cls.setNodeAttributes(node, path, int(level), dims, 0, 0, cls.lengthUnit(image), orientationSource, region)
        node.CreateDefaultDisplayNodes()
        attrs = (multiscales.root_attributes or {}).get("image-label") if multiscales.root_attributes else None
        if attrs is None:
            attrs = (readStoreAttributes(path) or {}).get("image-label")
        colorNode = cls.colorTableFromImageLabel(attrs, nodeName + "_colors")
        if colorNode is not None:
            node.GetDisplayNode().SetAndObserveColorNodeID(colorNode.GetID())
        cls.applyDisplayUnits(cls.lengthUnit(image))
        if Settings.get(Settings.LABELS_AS_SEGMENTATION, False):
            node = cls.segmentationFromLabelMap(node, colorNode)
        return [node]

    @staticmethod
    def segmentationFromLabelMap(labelNode, colorNode=None):
        """Replace a label map (and its colour table) by a Segmentation node carrying the same attributes."""
        segmentation = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", labelNode.GetName())
        segmentation.CreateDefaultDisplayNodes()
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(labelNode, segmentation)
        for key in labelNode.GetAttributeNames():
            segmentation.SetAttribute(key, labelNode.GetAttribute(key))
        slicer.mrmlScene.RemoveNode(labelNode)
        if colorNode is not None:
            slicer.mrmlScene.RemoveNode(colorNode)
        return segmentation

    @staticmethod
    def colorTableFromImageLabel(imageLabel, name):
        """Colour table node from ``image-label`` colours and properties, or None."""
        if not isinstance(imageLabel, dict):
            return None
        colors = {}
        for entry in imageLabel.get("colors") or []:
            try:
                value = int(entry["label-value"])
                rgba = [float(v) for v in entry.get("rgba", [0, 0, 0, 255])]
            except (KeyError, TypeError, ValueError):
                continue
            colors.setdefault(value, {})["rgba"] = rgba
        for entry in imageLabel.get("properties") or []:
            try:
                value = int(entry["label-value"])
            except (KeyError, TypeError, ValueError):
                continue
            label = entry.get("name") or entry.get("label") or entry.get("acronym")
            if label:
                colors.setdefault(value, {})["name"] = str(label)
        if not colors:
            return None
        maxValue = max(colors)
        if maxValue >= MAX_COLOR_TABLE_ENTRIES or min(colors) < 0:
            logging.warning(f"Label values up to {maxValue} are too many for a colour table; default colours are used.")
            return None
        colorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLColorTableNode", name)
        colorNode.SetTypeToUser()
        colorNode.SetHideFromEditors(False)
        colorNode.SetNumberOfColors(maxValue + 1)
        colorNode.SetColor(0, "background", 0.0, 0.0, 0.0, 0.0)
        for value, info in colors.items():
            r, g, b, a = (info.get("rgba") or [128.0, 128.0, 128.0, 255.0])[:4]
            colorNode.SetColor(value, info.get("name", str(value)), r / 255.0, g / 255.0, b / 255.0, a / 255.0)
        return colorNode

    # ---- region of interest ----

    @classmethod
    def regionFromRasBounds(cls, image, rasBounds):
        """Map RAS bounds [xmin,xmax,ymin,ymax,zmin,zmax] to index ranges at this level."""
        ijkToRas, _source = cls.ijkToRasMatrix(image)
        rasToIjk = np.linalg.inv(ijkToRas)
        corners = np.array(
            [[rasBounds[i], rasBounds[2 + j], rasBounds[4 + k], 1.0] for i in (0, 1) for j in (0, 1) for k in (0, 1)]
        )
        ijk = (rasToIjk @ corners.T).T[:, :3]
        shape = cls.spatialShape(image)
        region = {}
        for axisIndex, d in enumerate(SPATIAL_DIMS):
            if d not in shape:
                continue
            start = int(np.clip(np.floor(ijk[:, axisIndex].min() + 0.5), 0, shape[d]))
            stop = int(np.clip(np.ceil(ijk[:, axisIndex].max() + 0.5), 0, shape[d]))
            if stop <= start:
                raise ValueError("Region of interest does not intersect the image")
            region[d] = (start, stop)
        return region

    @classmethod
    def loadRegion(
        cls, path, roiNode, level=0, timeIndex=None, channels=None, name=None, userMessages=None, progress=None
    ):
        multiscales = cls.openMultiscales(path)
        bounds = [0.0] * 6
        roiNode.GetRASBounds(bounds)
        region = cls.regionFromRasBounds(multiscales.images[int(level)], bounds)
        return cls.loadImage(
            path,
            level=level,
            timeIndex=timeIndex,
            channels=channels,
            region=region,
            name=name,
            userMessages=userMessages,
            multiscales=multiscales,
            progress=progress,
        )

    # ---- view-driven refinement ----

    @staticmethod
    def sliceViewRasBounds(sliceViewName):
        """RAS bounds of the block a slice view shows: its field of view, extended along the
        normal by half the smaller side."""
        sliceWidget = slicer.app.layoutManager().sliceWidget(sliceViewName)
        if sliceWidget is None:
            raise ValueError(f"No slice view named '{sliceViewName}'")
        sliceNode = sliceWidget.mrmlSliceNode()
        sliceToRas = slicer.util.arrayFromVTKMatrix(sliceNode.GetSliceToRAS())
        width, height, _depth = sliceNode.GetFieldOfView()
        thickness = min(width, height) / 2.0
        corners = np.array(
            [
                [x, y, z, 1.0]
                for x in (-width / 2, width / 2)
                for y in (-height / 2, height / 2)
                for z in (-thickness / 2, thickness / 2)
            ]
        )
        ras = (sliceToRas @ corners.T).T[:, :3]
        return [ras[:, 0].min(), ras[:, 0].max(), ras[:, 1].min(), ras[:, 1].max(), ras[:, 2].min(), ras[:, 2].max()]

    @staticmethod
    def currentTimeIndex(path):
        """Selected item of the sequence browser showing this store, else 0."""
        for browser in slicer.util.getNodesByClass("vtkMRMLSequenceBrowserNode"):
            master = browser.GetMasterSequenceNode()
            if master is not None and samePath(master.GetAttribute("OMEZarr.Path"), path):
                return max(0, browser.GetSelectedItemNumber())
        return 0

    @staticmethod
    def sourceVolumeNode(path, channel, sliceViewName=None):
        """The full-extent volume loaded from this store and channel: the one shown as
        background of ``sliceViewName`` when it matches, else the most recently loaded."""

        def matches(node):
            return (
                node is not None
                and node.IsA("vtkMRMLScalarVolumeNode")
                and not node.IsA("vtkMRMLLabelMapVolumeNode")
                and samePath(node.GetAttribute("OMEZarr.Path"), path)
                and node.GetAttribute("OMEZarr.Channel") == str(channel)
                and not node.GetAttribute("OMEZarr.Refined")
                and not node.GetAttribute("OMEZarr.Region")
            )

        if sliceViewName:
            sliceWidget = slicer.app.layoutManager().sliceWidget(sliceViewName)
            if sliceWidget is not None:
                background = sliceWidget.sliceLogic().GetBackgroundLayer().GetVolumeNode()
                if matches(background):
                    return background
        candidates = [n for n in slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode") if matches(n)]
        return candidates[-1] if candidates else None

    @classmethod
    def matchDisplay(cls, node, path, sliceViewName=None):
        """Give a refined block the window/level and colour of the volume it refines."""
        source = cls.sourceVolumeNode(path, node.GetAttribute("OMEZarr.Channel"), sliceViewName)
        if source is None or source.GetDisplayNode() is None or node.GetDisplayNode() is None:
            return
        sourceDisplay, display = source.GetDisplayNode(), node.GetDisplayNode()
        display.SetAutoWindowLevel(False)
        display.SetWindowLevel(sourceDisplay.GetWindow(), sourceDisplay.GetLevel())
        if sourceDisplay.GetColorNodeID():
            display.SetAndObserveColorNodeID(sourceDisplay.GetColorNodeID())

    @classmethod
    def refinedNodes(cls, path=None, sliceViewName=None):
        nodes = slicer.util.getNodesByClass("vtkMRMLVolumeNode") + slicer.util.getNodesByClass(
            "vtkMRMLSegmentationNode"
        )
        return [
            n
            for n in nodes
            if n.GetAttribute("OMEZarr.Refined") == "1"
            and (path is None or samePath(n.GetAttribute("OMEZarr.Path"), path))
            and (sliceViewName is None or n.GetAttribute("OMEZarr.RefinedView") == sliceViewName)
        ]

    @classmethod
    def refineView(cls, path, sliceViewName="Red", maxBytes=None, timeIndex=None, userMessages=None, progress=None):
        """Reload what ``sliceViewName`` shows at the finest level whose block fits the budget.

        Each slice view keeps its own refined block, shown in that view's foreground layer
        (labels in its label layer) with the window/level of the coarse volume. The
        previous block of the same view and store is replaced. Returns the new nodes.
        """
        multiscales = cls.openMultiscales(path)
        bounds = cls.sliceViewRasBounds(sliceViewName)
        budget = maxBytes or cls.maxBytesFromSettings()
        channels = len(cls.channelDescriptions(multiscales, multiscales.images[0]))
        chosen = None
        for level, image in enumerate(multiscales.images):
            try:
                region = cls.regionFromRasBounds(image, bounds)
            except ValueError:
                continue
            if cls.volumeBytes(image, region) * channels <= budget:
                chosen = (level, region)
                break
        if chosen is None:
            raise ValueError("The view does not intersect the image, or no level fits the memory budget")
        level, region = chosen
        shown = cls.sourceVolumeNode(path, 0, sliceViewName)
        shownLevel = shown.GetAttribute("OMEZarr.Level") if shown is not None else None
        if shownLevel is not None and level >= int(shownLevel):
            raise ValueError("Zoom in: no finer level than the one displayed fits the memory budget for this view")
        for node in cls.refinedNodes(path, sliceViewName):
            display = node.GetDisplayNode() if node.IsA("vtkMRMLVolumeNode") else None
            colorNode = display.GetColorNode() if display else None
            slicer.mrmlScene.RemoveNode(node)
            if colorNode is not None and colorNode.GetAttribute("OMEZarr.Refined"):
                slicer.mrmlScene.RemoveNode(colorNode)
        nodes = cls.loadImage(
            path,
            level=level,
            timeIndex=cls.currentTimeIndex(path) if timeIndex is None else timeIndex,
            region=region,
            name=f"{cls.defaultNodeName(path)}_{sliceViewName}",
            userMessages=userMessages,
            multiscales=multiscales,
            progress=progress,
        )
        scalars, labels = [], []
        for node in nodes:
            node.SetAttribute("OMEZarr.Refined", "1")
            node.SetAttribute("OMEZarr.RefinedView", sliceViewName)
            if node.IsA("vtkMRMLLabelMapVolumeNode"):
                labels.append(node)
                colorNode = node.GetDisplayNode().GetColorNode() if node.GetDisplayNode() else None
                if colorNode is not None and colorNode.GetAttribute("OMEZarr.Path") is None:
                    colorNode.SetAttribute("OMEZarr.Refined", "1")
            elif node.IsA("vtkMRMLScalarVolumeNode"):
                scalars.append(node)
                cls.matchDisplay(node, path, sliceViewName)
        composite = slicer.app.layoutManager().sliceWidget(sliceViewName).sliceLogic().GetSliceCompositeNode()
        if scalars:
            composite.SetForegroundVolumeID(scalars[0].GetID())
            composite.SetForegroundOpacity(1.0)
        if labels:
            composite.SetLabelVolumeID(labels[0].GetID())
        return nodes

    # ---- automatic refinement ----

    _autoRefiners = {}

    DEFAULT_AUTO_REFINE_VIEWS = ("Red", "Yellow", "Green")

    @classmethod
    def startAutoRefine(cls, path, sliceViewNames=DEFAULT_AUTO_REFINE_VIEWS, delayMs=600, maxBytes=None):
        """Refine each of ``sliceViewNames`` whenever it stops moving. Returns the AutoRefiner.

        The memory budget is shared equally between the views.
        """
        cls.stopAutoRefine(path)
        if isinstance(sliceViewNames, str):
            sliceViewNames = (sliceViewNames,)
        refiner = AutoRefiner(str(path), list(sliceViewNames), delayMs, maxBytes)
        cls._autoRefiners[str(path)] = refiner
        return refiner

    @classmethod
    def stopAutoRefine(cls, path=None):
        keys = [str(path)] if path is not None else list(cls._autoRefiners)
        for key in keys:
            refiner = cls._autoRefiners.pop(key, None)
            if refiner is not None:
                refiner.stop()

    @classmethod
    def autoRefiner(cls, path):
        return cls._autoRefiners.get(str(path))

    # ---- display units ----

    @staticmethod
    def applyDisplayUnits(lengthUnit):
        """Switch Slicer's length display to the store's unit when the setting is on."""
        symbol = UNIT_SYMBOLS.get(lengthUnit)
        if not Settings.get(Settings.DISPLAY_UNITS, False) or symbol in (None, "mm"):
            return
        unitNode = slicer.mrmlScene.GetNodeByID("vtkMRMLUnitNodeApplicationLength")
        if unitNode is None:
            return
        unitNode.SetDisplayCoefficient(1.0 / LENGTH_UNIT_TO_MM[lengthUnit])
        unitNode.SetSuffix(symbol)
        unitNode.SetPrecision(2)

    @staticmethod
    def resetDisplayUnits():
        unitNode = slicer.mrmlScene.GetNodeByID("vtkMRMLUnitNodeApplicationLength")
        if unitNode is None:
            return
        unitNode.SetDisplayCoefficient(1.0)
        unitNode.SetSuffix("mm")
        unitNode.SetPrecision(3)

    # ---- write side ----

    @staticmethod
    def ngffImageFromVolumeNode(volumeNode, name=None):
        """NgffImage (z,y,x, millimeter, RFC-4 orientation) from a scalar or label map volume node."""
        import ngff_zarr
        from ngff_zarr.rfc4 import itk_direction_to_anatomical_orientation

        array = slicer.util.arrayFromVolume(volumeNode)
        if array.ndim != 3:
            raise ValueError("Only single-component volumes can be written as OME-Zarr")
        ijkToRasVtk = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijkToRasVtk)
        ijkToRas = slicer.util.arrayFromVTKMatrix(ijkToRasVtk)
        ijkToLps = np.diag([-1.0, -1.0, 1.0, 1.0]) @ ijkToRas
        spacing = np.linalg.norm(ijkToLps[:3, :3], axis=0)
        direction = ijkToLps[:3, :3] / spacing
        origin = ijkToLps[:3, 3]
        orientations = {
            d: itk_direction_to_anatomical_orientation(list(direction[:, i])) for i, d in enumerate(SPATIAL_DIMS)
        }
        image = ngff_zarr.to_ngff_image(
            np.ascontiguousarray(array),
            dims=("z", "y", "x"),
            scale={d: float(spacing[i]) for i, d in enumerate(SPATIAL_DIMS)},
            translation={d: float(origin[i]) for i, d in enumerate(SPATIAL_DIMS)},
            name=name or volumeNode.GetName(),
            axes_units={d: "millimeter" for d in SPATIAL_DIMS},
        )
        image.axes_orientations = orientations
        return image

    @staticmethod
    def imageLabelFromNode(labelNode, version):
        """``image-label`` metadata for a label map volume, from its colour table."""
        values = np.unique(slicer.util.arrayFromVolume(labelNode))
        values = [int(v) for v in values if v != 0]
        colors, properties = [], []
        display = labelNode.GetDisplayNode()
        colorNode = display.GetColorNode() if display else None
        for value in values:
            rgba = [0.5, 0.5, 0.5, 1.0]
            name = str(value)
            if colorNode is not None and value < colorNode.GetNumberOfColors():
                colorNode.GetColor(value, rgba)
                name = colorNode.GetColorName(value) or name
            colors.append({"label-value": value, "rgba": [int(round(c * 255)) for c in rgba]})
            properties.append({"label-value": value, "name": name})
        return {"version": version, "colors": colors, "properties": properties}

    @classmethod
    def writeVolume(cls, node, storePath, progress=None):
        """Write a scalar volume, label map or segmentation as an OME-Zarr multiscales store."""
        ngff_zarr = cls.ensureNgffZarr()
        from ngff_zarr import Methods

        if isRemoteUrl(storePath):
            raise ValueError("Writing to a remote store is not supported; write to a local directory and upload it.")
        if node.IsA("vtkMRMLSegmentationNode"):
            labelNode = cls.labelMapFromSegmentation(node)
            try:
                return cls.writeVolume(labelNode, storePath, progress)
            finally:
                display = labelNode.GetDisplayNode()
                colorNode = display.GetColorNode() if display else None
                slicer.mrmlScene.RemoveNode(labelNode)
                if colorNode is not None:
                    slicer.mrmlScene.RemoveNode(colorNode)
        isLabel = node.IsA("vtkMRMLLabelMapVolumeNode")
        image = cls.ngffImageFromVolumeNode(node)
        method = Methods.ITKWASM_LABEL_IMAGE if isLabel else None

        def write():
            try:
                multiscales = ngff_zarr.to_multiscales(image, method=method)
            except Exception:  # noqa: BLE001 - fall back to a single level rather than fail
                logging.exception("Multiscale generation failed; writing a single level")
                multiscales = ngff_zarr.to_multiscales(image, scale_factors=[])
            ngff_zarr.to_ome_zarr(storePath, multiscales, overwrite=True)
            return multiscales

        if progress:
            progress(0, 1, f"Writing {os.path.basename(storePath)}")
        multiscales = runResponsive(write)
        version = str(getattr(multiscales.metadata, "version", "0.5") or "0.5")
        if isLabel:
            cls.addImageLabelMetadata(storePath, cls.imageLabelFromNode(node, version))
            cls.registerLabelInParent(storePath)
        cls.clearCache()
        if progress:
            progress(1, 1)
        return storePath

    @staticmethod
    def labelMapFromSegmentation(segmentationNode):
        """Temporary label map (with a colour table of segment names and colours) of all segments."""
        labelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", segmentationNode.GetName())
        logic = slicer.modules.segmentations.logic()
        if not logic.ExportAllSegmentsToLabelmapNode(
            segmentationNode, labelNode, slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY
        ):
            slicer.mrmlScene.RemoveNode(labelNode)
            raise ValueError("Could not export the segmentation to a label map")
        return labelNode

    @staticmethod
    def _updateGroupMetadata(groupPath, update):
        """Apply ``update(omeAttrs)`` to the OME attributes of a local Zarr group (v2 or v3)."""
        v3 = os.path.join(groupPath, "zarr.json")
        v2 = os.path.join(groupPath, ".zattrs")
        if os.path.isfile(v3):
            with open(v3, encoding="utf-8") as fp:
                document = json.load(fp)
            attributes = document.setdefault("attributes", {})
            ome = attributes.setdefault("ome", {})
            update(ome)
            with open(v3, "w", encoding="utf-8") as fp:
                json.dump(document, fp, indent=2)
        else:
            document = {}
            if os.path.isfile(v2):
                with open(v2, encoding="utf-8") as fp:
                    document = json.load(fp)
            update(document)
            with open(v2, "w", encoding="utf-8") as fp:
                json.dump(document, fp, indent=2)

    @classmethod
    def addImageLabelMetadata(cls, storePath, imageLabel):
        cls._updateGroupMetadata(storePath, lambda ome: ome.__setitem__("image-label", imageLabel))

    @classmethod
    def registerLabelInParent(cls, storePath):
        """When written into ``<image>/labels/<name>``, list the name in the ``labels`` group."""
        labelsDir = os.path.dirname(os.path.abspath(storePath))
        imageRoot = os.path.dirname(labelsDir)
        if os.path.basename(labelsDir) != "labels" or not omeZarrRootFromPath(imageRoot):
            return
        name = os.path.basename(storePath)
        zarrV3 = os.path.isfile(os.path.join(imageRoot, "zarr.json"))
        version = (readStoreAttributes(imageRoot) or {}).get("version", "0.4")
        if zarrV3 and not os.path.isfile(os.path.join(labelsDir, "zarr.json")):
            with open(os.path.join(labelsDir, "zarr.json"), "w", encoding="utf-8") as fp:
                json.dump({"zarr_format": 3, "node_type": "group", "attributes": {"ome": {"version": version}}}, fp)
        elif not zarrV3 and not os.path.isfile(os.path.join(labelsDir, ".zgroup")):
            with open(os.path.join(labelsDir, ".zgroup"), "w", encoding="utf-8") as fp:
                json.dump({"zarr_format": 2}, fp)

        def update(ome):
            names = [n for n in ome.get("labels", []) if n != name]
            ome["labels"] = [*names, name]

        cls._updateGroupMetadata(labelsDir, update)


class AutoRefiner:
    """Reloads the block each observed slice view shows once that view has been still.

    Every change of a slice node restarts a single timer, so nothing loads while the user
    pans or zooms. On timeout, only the views whose bounds changed are refined.
    """

    def __init__(self, path, sliceViewNames, delayMs, maxBytes):
        self.path = path
        self.sliceViewNames = list(sliceViewNames)
        self.maxBytes = maxBytes
        self.lastBounds = {}
        self.busy = False
        self.pending = False
        self.refreshCount = 0
        self.timer = qt.QTimer()
        self.timer.setSingleShot(True)
        self.timer.setInterval(int(delayMs))
        self.timer.timeout.connect(self.refresh)
        self.observers = []
        layoutManager = slicer.app.layoutManager()
        for viewName in self.sliceViewNames:
            sliceWidget = layoutManager.sliceWidget(viewName)
            if sliceWidget is None:
                raise ValueError(f"No slice view named '{viewName}'")
            sliceNode = sliceWidget.mrmlSliceNode()
            self.observers.append(
                (sliceNode, sliceNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.onSliceModified))
            )
        self.sceneObserverTag = slicer.mrmlScene.AddObserver(slicer.vtkMRMLScene.EndCloseEvent, self.onSceneClosed)
        self.timer.start()

    def stop(self):
        self.timer.stop()
        for node, tag in self.observers:
            node.RemoveObserver(tag)
        self.observers = []
        if self.sceneObserverTag is not None:
            slicer.mrmlScene.RemoveObserver(self.sceneObserverTag)
            self.sceneObserverTag = None

    def onSceneClosed(self, caller=None, event=None):
        OMEZarrLogic.stopAutoRefine(self.path)

    def onSliceModified(self, caller=None, event=None):
        if self.busy:
            self.pending = True  # the application stays responsive while a block loads
        else:
            self.timer.start()

    @staticmethod
    def tolerance(bounds):
        """A view counts as still when it moved by less than 1% of its smallest extent."""
        return 0.01 * min(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])

    def idle(self):
        return not self.busy and not self.timer.isActive()

    def viewBudget(self):
        return (self.maxBytes or OMEZarrLogic.maxBytesFromSettings()) // max(1, len(self.sliceViewNames))

    def refresh(self):
        if self.busy:
            self.timer.start()
            return
        self.busy = True
        try:
            for viewName in self.sliceViewNames:
                try:
                    bounds = OMEZarrLogic.sliceViewRasBounds(viewName)
                except ValueError:
                    continue
                previous = self.lastBounds.get(viewName)
                if previous is not None and np.allclose(bounds, previous, rtol=0.0, atol=self.tolerance(bounds)):
                    continue
                self.lastBounds[viewName] = bounds
                try:
                    OMEZarrLogic.refineView(self.path, viewName, maxBytes=self.viewBudget())
                    self.refreshCount += 1
                except (ValueError, InterruptedError) as e:
                    logging.debug(f"Auto-refine of {viewName} skipped: {e}")
                except Exception:  # noqa: BLE001 - never let a timer callback raise into Qt
                    logging.exception(f"Auto-refine of {viewName} failed")
        finally:
            self.busy = False
            if self.pending:
                self.pending = False
                self.timer.start()


#
# File reader (drives Add Data, drag-and-drop, slicer.util.loadNodeFromFile)
#


class OMEZarrFileReader:
    def __init__(self, parent):
        self.parent = parent

    def description(self):
        return _("OME-Zarr image")

    def fileType(self):
        return "OMEZarr"

    def extensions(self):
        return [
            _("OME-Zarr") + " (*.ome.zarr)",
            _("OME-Zarr") + " (*.zarr)",
            _("OME-Zarr zip") + " (*.ozx *.zarr.zip)",
            _("OME-Zarr metadata") + " (zarr.json .zattrs)",
        ]

    def canLoadFileConfidence(self, filePath):
        # Do not use self.parent.supportedNameFilters(): it rejects directories.
        # 0.9 outranks the default extension-based confidence of core readers
        # (0.5 + 0.01 * extension length), which matters for zarr.json/.zattrs
        # entries listed by the Add Data dialog.
        return 0.9 if omeZarrRootFromPath(filePath) else 0.0

    def load(self, properties):
        try:
            root = omeZarrRootFromPath(properties["fileName"])
            if not root:
                raise ValueError(f"Not an OME-Zarr multiscales store: {properties['fileName']}")

            def optional(key, cast):
                value = properties.get(key)
                return cast(value) if value not in (None, "") else None

            with Progress(_("Loading OME-Zarr...")) as progress:
                nodes = OMEZarrLogic.loadImage(
                    root,
                    level=optional("level", int),
                    timeIndex=optional("timeIndex", int),
                    name=properties.get("name"),
                    maxBytes=optional("maxBytes", int),
                    labels=optional("labels", lambda v: str(v).lower() in ("true", "1")),
                    timeMode=optional("timeMode", str),
                    asLabelMap=optional("asLabelMap", lambda v: str(v).lower() in ("true", "1")),
                    userMessages=self.parent.userMessages(),
                    progress=progress,
                )
        except InterruptedError:
            self.parent.userMessages().AddMessage(vtk.vtkCommand.WarningEvent, "Loading cancelled")
            return False
        except Exception as e:  # noqa: BLE001 - report everything to the user
            import traceback

            traceback.print_exc()
            self.parent.userMessages().AddMessage(vtk.vtkCommand.ErrorEvent, f"Failed to read OME-Zarr: {e}")
            return False

        if properties.get("show", True) and nodes:
            scalars = scalarVolumes(nodes)
            labels = labelMaps(nodes)
            slicer.util.setSliceViewerLayers(
                background=scalars[0] if scalars else "keep-current",
                label=labels[0] if labels else "keep-current",
                fit=True,
            )
            coarse = any(n.GetAttribute("OMEZarr.Level") not in (None, "0") for n in scalars)
            if coarse and Settings.get(Settings.AUTO_REFINE, False) and slicer.util.mainWindow():
                OMEZarrLogic.startAutoRefine(root)
        self.parent.loadedNodes = [node.GetID() for node in nodes]
        return True


#
# File writer (File > Save, slicer.util.saveNode)
#


class OMEZarrFileWriter:
    def __init__(self, parent):
        self.parent = parent

    def description(self):
        return _("OME-Zarr image")

    def fileType(self):
        return "OMEZarr"

    def extensions(self, obj):
        # The second, generic filter lets a label map be written as ``<image>/labels/<name>``.
        return [_("OME-Zarr") + " (.ome.zarr)", _("OME-Zarr directory") + " (*)"]

    def canWriteObjectConfidence(self, obj):
        # Below the default 0.5 so the native formats stay the default choice in the save dialog.
        if (
            obj.IsA("vtkMRMLLabelMapVolumeNode")
            or obj.IsA("vtkMRMLScalarVolumeNode")
            or obj.IsA("vtkMRMLSegmentationNode")
        ):
            return 0.3
        return 0.0

    def write(self, properties):
        try:
            node = slicer.mrmlScene.GetNodeByID(properties["nodeID"])
            with Progress(_("Writing OME-Zarr...")) as progress:
                OMEZarrLogic.writeVolume(node, properties["fileName"], progress)
        except Exception as e:  # noqa: BLE001 - report everything to the user
            import traceback

            traceback.print_exc()
            self.parent.userMessages().AddMessage(vtk.vtkCommand.ErrorEvent, f"Failed to write OME-Zarr: {e}")
            return False
        self.parent.writtenNodes = [node.GetID()]
        return True


#
# Drop target for OME-Zarr directories (same mechanism the DICOM module uses)
#


class OMEZarrFileDialog:
    """Detected by name by qSlicerScriptedLoadableModule and registered with the IO manager."""

    def __init__(self, qSlicerFileDialog):
        self.qSlicerFileDialog = qSlicerFileDialog
        qSlicerFileDialog.fileType = "OMEZarr"
        qSlicerFileDialog.description = _("Load OME-Zarr image")
        qSlicerFileDialog.action = slicer.qSlicerFileDialog.Read
        self.pathsToLoad = []

    def execDialog(self):
        # Not used: loading is triggered from dropEvent.
        return True

    def isMimeDataAccepted(self):
        """Accept only when every dropped URL is an OME-Zarr store (so plain files/dirs keep their usual handling)."""
        self.pathsToLoad = []
        mimeData = self.qSlicerFileDialog.mimeData()
        accepted = False
        if mimeData.hasFormat("text/uri-list"):
            roots = [omeZarrRootFromPath(url.toLocalFile() or url.toString()) for url in mimeData.urls()]
            if roots and all(roots):
                self.pathsToLoad = roots
                accepted = True
        self.qSlicerFileDialog.acceptMimeData(accepted)

    def dropEvent(self):
        paths = list(self.pathsToLoad)
        self.pathsToLoad = []
        # Defer so the drag source application is not blocked while we read data.
        qt.QTimer.singleShot(0, lambda: self._loadPaths(paths))

    @staticmethod
    def _loadPaths(paths):
        for path in paths:
            with slicer.util.tryWithErrorDisplay(_("Failed to load OME-Zarr image"), waitCursor=True):
                slicer.util.loadNodeFromFile(path, "OMEZarr")


#
# Widget
#


class OMEZarrWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        import ctk

        self.logic = OMEZarrLogic()
        self.multiscales = None
        self.path = None
        self._inspecting = False

        # -- Store --
        storeBox = ctk.ctkCollapsibleButton()
        storeBox.text = _("Store")
        self.layout.addWidget(storeBox)
        storeLayout = qt.QVBoxLayout(storeBox)  # a form layout mis-measures word-wrapped labels

        self.pathEdit = ctk.ctkPathLineEdit()
        self.pathEdit.filters = ctk.ctkPathLineEdit.Dirs
        self.pathEdit.settingKey = "OMEZarr/LastPath"
        self.pathEdit.setToolTip(_("Local .ome.zarr directory, .ozx file, or https:// / s3:// URL"))
        self.inspectButton = qt.QToolButton()
        self.inspectButton.setIcon(slicer.app.style().standardIcon(qt.QStyle.SP_BrowserReload))
        self.inspectButton.setToolTip(_("Read the store's metadata again"))
        pathRow = qt.QHBoxLayout()
        pathRow.addWidget(self.pathEdit, 1)
        pathRow.addWidget(self.inspectButton)
        storeLayout.addLayout(pathRow)

        self.infoLabel = qt.QLabel(
            _("Choose an OME-Zarr store, or drop one onto Slicer: its resolution levels are listed here.")
        )
        self.infoLabel.wordWrap = True
        self.infoLabel.setTextInteractionFlags(qt.Qt.TextSelectableByMouse)
        storeLayout.addWidget(self.infoLabel)

        self.levelTable = qt.QTableWidget(0, 4)
        self.levelTable.setHorizontalHeaderLabels([_("Level"), _("Voxels (x, y, z)"), _("Spacing"), _("Memory")])
        self.levelTable.setTextElideMode(qt.Qt.ElideRight)
        self.levelTable.setWordWrap(False)
        self.levelTable.verticalHeader().setVisible(False)
        self.levelTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.levelTable.setSelectionMode(qt.QAbstractItemView.SingleSelection)
        self.levelTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.levelTable.setHorizontalScrollBarPolicy(qt.Qt.ScrollBarAlwaysOff)
        self.levelTable.setVerticalScrollBarPolicy(qt.Qt.ScrollBarAlwaysOff)
        header = self.levelTable.horizontalHeader()
        header.setSectionResizeMode(qt.QHeaderView.Stretch)
        for column in (0, 1, 3):
            header.setSectionResizeMode(column, qt.QHeaderView.ResizeToContents)
        storeLayout.addWidget(self.levelTable)

        self.legendLabel = qt.QLabel(_("Bold: the level the memory budget selects. ✓: loaded in the scene."))
        self.legendLabel.wordWrap = True
        self.legendLabel.enabled = False  # greyed like a hint, in any Slicer style
        storeLayout.addWidget(self.legendLabel)

        self.timeIndexLabel = qt.QLabel(_("Time point:"))
        self.timeIndexSpinBox = qt.QSpinBox()
        self.timeIndexSpinBox.setRange(0, 0)
        self.timeIndexSpinBox.setSpecialValueText(_("all (sequence)"))
        timeRow = qt.QHBoxLayout()
        timeRow.addWidget(self.timeIndexLabel)
        timeRow.addWidget(self.timeIndexSpinBox, 1)
        storeLayout.addLayout(timeRow)

        self.loadButton = qt.QPushButton(_("Load selected level"))
        storeLayout.addWidget(self.loadButton)

        # -- Full resolution --
        refineBox = ctk.ctkCollapsibleButton()
        refineBox.text = _("Full resolution")
        self.layout.addWidget(refineBox)
        refineBoxLayout = qt.QVBoxLayout(refineBox)
        refineLayout = qt.QFormLayout()
        refineBoxLayout.addLayout(refineLayout)

        self.viewSelector = qt.QComboBox()
        self.viewSelector.setToolTip(_("The 2D view that 'Refine view' and 'New ROI in view' act on"))
        self.updateViewSelector()
        self.refineButton = qt.QPushButton(_("Refine view"))
        self.refineButton.setToolTip(
            _("Reload the block shown by the slice view at the finest level that fits the memory budget")
        )
        refineRow = qt.QHBoxLayout()
        refineRow.addWidget(self.viewSelector)
        refineRow.addWidget(self.refineButton, 1)
        refineLayout.addRow(_("2D view:"), refineRow)

        self.autoRefineCheckBox = qt.QCheckBox(_("Refine the slice views automatically while browsing"))
        self.autoRefineCheckBox.setToolTip(_("Reloads a view's block after it has been still for half a second"))
        refineLayout.addRow(self.autoRefineCheckBox)

        self.roiSelector = slicer.qMRMLNodeComboBox()
        self.roiSelector.nodeTypes = ["vtkMRMLMarkupsROINode"]
        self.roiSelector.addEnabled = False  # "New ROI in view" creates one that is already placed
        self.roiSelector.removeEnabled = True
        self.roiSelector.renameEnabled = True
        self.roiSelector.noneEnabled = True
        self.roiSelector.setMRMLScene(slicer.mrmlScene)
        self.createRoiButton = qt.QPushButton(_("New ROI in view"))
        self.createRoiButton.setToolTip(
            _("Create a region of interest covering the middle of the selected slice view; drag its handles to adjust")
        )
        roiRow = qt.QHBoxLayout()
        roiRow.addWidget(self.roiSelector, 1)
        roiRow.addWidget(self.createRoiButton)
        refineLayout.addRow(_("Region of interest:"), roiRow)
        self.loadRegionButton = qt.QPushButton(_("Load region at the selected level"))
        self.loadRegionButton.setToolTip(_("Only the chunks the region intersects are read"))
        refineLayout.addRow(self.loadRegionButton)

        self.statusLabel = qt.QLabel()
        self.statusLabel.wordWrap = True
        refineBoxLayout.addWidget(self.statusLabel)

        # -- Settings --
        settingsBox = ctk.ctkCollapsibleButton()
        settingsBox.text = _("Settings")
        settingsBox.collapsed = True
        self.layout.addWidget(settingsBox)
        settingsLayout = qt.QFormLayout(settingsBox)

        self.maxBytesSpinBox = qt.QSpinBox()
        self.maxBytesSpinBox.setRange(0, 1 << 20)
        self.maxBytesSpinBox.setSuffix(" MiB")
        self.maxBytesSpinBox.setSpecialValueText(_("automatic (25% of free RAM)"))
        self.maxBytesSpinBox.setValue(Settings.get(Settings.MAX_BYTES, 0) >> 20)
        self.maxBytesSpinBox.setToolTip(_("Budget for automatic level selection"))
        settingsLayout.addRow(_("Memory budget:"), self.maxBytesSpinBox)

        self.orientationSelector = qt.QComboBox()
        self.orientationSelector.addItems(["LPS", "RAS"])
        self.orientationSelector.setCurrentText(Settings.get(Settings.ORIENTATION, "LPS"))
        self.orientationSelector.setToolTip(_("How x/y/z axes are interpreted when a store has no RFC-4 orientation"))
        settingsLayout.addRow(_("Axes without RFC-4:"), self.orientationSelector)

        self.loadLabelsCheckBox = qt.QCheckBox()
        self.loadLabelsCheckBox.checked = Settings.get(Settings.LOAD_LABELS, True)
        settingsLayout.addRow(_("Load labels:"), self.loadLabelsCheckBox)

        self.labelsAsSegmentationCheckBox = qt.QCheckBox()
        self.labelsAsSegmentationCheckBox.checked = Settings.get(Settings.LABELS_AS_SEGMENTATION, False)
        self.labelsAsSegmentationCheckBox.setToolTip(_("Import labels as Segmentation nodes instead of label maps"))
        settingsLayout.addRow(_("Labels as segmentation:"), self.labelsAsSegmentationCheckBox)

        self.timeModeSelector = qt.QComboBox()
        self.timeModeSelector.addItem(_("all time points as a Sequence"), "sequence")
        self.timeModeSelector.addItem(_("first time point only"), "index")
        self.timeModeSelector.setCurrentIndex(0 if Settings.get(Settings.TIME_MODE, "sequence") == "sequence" else 1)
        settingsLayout.addRow(_("Time axis:"), self.timeModeSelector)

        self.displayUnitsCheckBox = qt.QCheckBox()
        self.displayUnitsCheckBox.checked = Settings.get(Settings.DISPLAY_UNITS, False)
        self.displayUnitsCheckBox.setToolTip(_("Show lengths in the store's unit (µm, nm) instead of mm"))
        settingsLayout.addRow(_("Display in store units:"), self.displayUnitsCheckBox)

        self.detectLabelMapsCheckBox = qt.QCheckBox()
        self.detectLabelMapsCheckBox.checked = Settings.get(Settings.DETECT_LABEL_MAPS, True)
        self.detectLabelMapsCheckBox.setToolTip(_("Integer stores with few distinct values load as label maps"))
        settingsLayout.addRow(_("Detect label maps:"), self.detectLabelMapsCheckBox)

        self.storageOptionsEdit = qt.QLineEdit(Settings.get(Settings.STORAGE_OPTIONS, ""))
        self.storageOptionsEdit.setPlaceholderText('{"anon": true, "region": "us-west-2"}')
        self.storageOptionsEdit.setToolTip(
            _("JSON storage options for remote stores (S3 credentials, endpoint, region)")
        )
        settingsLayout.addRow(_("Remote storage options:"), self.storageOptionsEdit)

        self.autoRefineOnLoadCheckBox = qt.QCheckBox()
        self.autoRefineOnLoadCheckBox.checked = Settings.get(Settings.AUTO_REFINE, False)
        self.autoRefineOnLoadCheckBox.setToolTip(
            _("After loading a downsampled level, refine the slice views automatically")
        )
        settingsLayout.addRow(_("Auto-refine after load:"), self.autoRefineOnLoadCheckBox)

        self.resetUnitsButton = qt.QPushButton(_("Reset display to mm"))
        settingsLayout.addRow(self.resetUnitsButton)

        self.layout.addStretch(1)

        self.inspectButton.connect("clicked(bool)", self.onInspect)
        self.pathEdit.connect("currentPathChanged(QString)", self.onPathChanged)
        self.loadButton.connect("clicked(bool)", self.onLoad)
        self.refineButton.connect("clicked(bool)", self.onRefine)
        self.autoRefineCheckBox.connect("toggled(bool)", self.onAutoRefineToggled)
        self.loadRegionButton.connect("clicked(bool)", self.onLoadRegion)
        self.createRoiButton.connect("clicked(bool)", self.onCreateRoi)
        self.levelTable.connect("itemSelectionChanged()", self.updateButtons)
        self.levelTable.connect("cellDoubleClicked(int,int)", lambda row, column: self.onLoad())
        self.roiSelector.connect("currentNodeChanged(vtkMRMLNode*)", lambda node: self.updateButtons())
        self.updateButtons()
        self.maxBytesSpinBox.connect("valueChanged(int)", lambda mib: Settings.set(Settings.MAX_BYTES, int(mib) << 20))
        self.orientationSelector.connect("currentTextChanged(QString)", lambda t: Settings.set(Settings.ORIENTATION, t))
        self.loadLabelsCheckBox.connect("toggled(bool)", lambda b: Settings.set(Settings.LOAD_LABELS, bool(b)))
        self.labelsAsSegmentationCheckBox.connect(
            "toggled(bool)", lambda b: Settings.set(Settings.LABELS_AS_SEGMENTATION, bool(b))
        )
        self.timeModeSelector.connect(
            "currentIndexChanged(int)", lambda i: Settings.set(Settings.TIME_MODE, self.timeModeSelector.itemData(i))
        )
        self.displayUnitsCheckBox.connect("toggled(bool)", lambda b: Settings.set(Settings.DISPLAY_UNITS, bool(b)))
        self.resetUnitsButton.connect("clicked(bool)", OMEZarrLogic.resetDisplayUnits)
        self.autoRefineOnLoadCheckBox.connect("toggled(bool)", lambda b: Settings.set(Settings.AUTO_REFINE, bool(b)))
        self.detectLabelMapsCheckBox.connect(
            "toggled(bool)", lambda b: Settings.set(Settings.DETECT_LABEL_MAPS, bool(b))
        )
        self.storageOptionsEdit.connect(
            "editingFinished()", lambda: Settings.set(Settings.STORAGE_OPTIONS, self.storageOptionsEdit.text)
        )

    def updateViewSelector(self):
        """List the slice views of the current layout as "Red (Axial)", keeping the selection."""
        layoutManager = slicer.app.layoutManager()
        current = self.viewSelector.currentData or "Red"
        self.viewSelector.clear()
        for name in layoutManager.sliceViewNames():
            orientation = layoutManager.sliceWidget(name).mrmlSliceNode().GetOrientation()
            self.viewSelector.addItem(f"{name} ({orientation})", name)
        self.viewSelector.setCurrentIndex(max(0, self.viewSelector.findData(current)))

    def enter(self):
        """Show the store of the displayed OME-Zarr volume, else the last one used."""
        self.updateViewSelector()
        if self.path:
            self.updateLevelStatus()
            return
        sliceWidget = slicer.app.layoutManager().sliceWidget("Red")
        background = sliceWidget.sliceLogic().GetBackgroundLayer().GetVolumeNode() if sliceWidget else None
        candidates = [background, *slicer.util.getNodesByClass("vtkMRMLVolumeNode")[::-1]]
        stored = next(
            (n.GetAttribute("OMEZarr.Path") for n in candidates if n and n.GetAttribute("OMEZarr.Path")), None
        )
        if stored and not samePath(stored, self.pathEdit.currentPath):
            self.pathEdit.currentPath = stored  # triggers onPathChanged
        else:
            self.onPathChanged(self.pathEdit.currentPath)

    def onPathChanged(self, path):
        """Inspect as soon as the path names a store; stay quiet while it does not."""
        if self._inspecting:
            return
        if omeZarrRootFromPath(path):
            if not samePath(path, self.path):
                self.onInspect()
        else:
            self.path = None
            self.multiscales = None
            self.levelTable.setRowCount(0)
            self.fitTableHeight()
            self.infoLabel.text = _(
                "Choose an OME-Zarr store, or drop one onto Slicer: its resolution levels are listed here."
            )
            self.updateButtons()

    def selectedLevel(self):
        rows = self.levelTable.selectionModel().selectedRows() if self.levelTable.rowCount else []
        return rows[0].row() if rows else -1

    def updateButtons(self):
        hasStore = self.path is not None
        self.loadButton.setEnabled(hasStore and self.selectedLevel() >= 0)
        self.refineButton.setEnabled(hasStore)
        self.createRoiButton.setEnabled(hasStore)
        self.autoRefineCheckBox.setEnabled(hasStore)
        self.loadRegionButton.setEnabled(hasStore and self.roiSelector.currentNode() is not None)

    def fitTableHeight(self):
        """Size the table to its rows so that it never shows an empty area."""
        table = self.levelTable
        rows = sum(table.rowHeight(row) for row in range(table.rowCount)) or table.verticalHeader().defaultSectionSize
        table.setFixedHeight(table.horizontalHeader().height + rows + 2 * table.frameWidth)

    def onInspect(self):
        path = omeZarrRootFromPath(self.pathEdit.currentPath)
        if not path:
            slicer.util.errorDisplay(_("Not an OME-Zarr multiscales store."))
            return
        self._inspecting = True  # adding the path to the history re-emits currentPathChanged
        try:
            with slicer.util.tryWithErrorDisplay(_("Failed to open store"), waitCursor=True):
                self.multiscales = self.logic.openMultiscales(path, useCache=False)
                self.path = path
                self.pathEdit.addCurrentPathToHistory()
                self.showLevels()
        finally:
            self._inspecting = False
        self.updateButtons()

    def showLevels(self):
        """Fill the level table and the summary line from the open store."""
        info = self.logic.levelInfo(self.multiscales)
        self.levelTable.setRowCount(len(info))
        for row, level in enumerate(info):
            shape = dict(zip(level["dims"], level["shape"], strict=False))
            unit = next((level["units"].get(d) for d in SPATIAL_DIMS if level["units"].get(d)), None)
            unit = UNIT_SYMBOLS.get(unit, unit or "")
            values = [
                str(level["level"]),
                " × ".join(str(shape[d]) for d in SPATIAL_DIMS if d in shape),
                (" × ".join(f"{level['scale'][d]:.4g}" for d in SPATIAL_DIMS if d in shape) + f" {unit}").strip(),
                self.formatBytes(level["bytes"]),
            ]
            for column, value in enumerate(values):
                item = qt.QTableWidgetItem(value)
                item.setToolTip(value)
                if column == 3:
                    item.setTextAlignment(qt.Qt.AlignRight | qt.Qt.AlignVCenter)
                self.levelTable.setItem(row, column, item)
        self.fitTableHeight()

        image = self.multiscales.images[0]
        channels = self.logic.axisLength(image, "c")
        timePoints = self.logic.axisLength(image, "t")
        self.timeIndexSpinBox.setRange(-1 if timePoints > 1 else 0, max(0, timePoints - 1))
        self.timeIndexSpinBox.value = -1 if timePoints > 1 else 0
        self.timeIndexLabel.setVisible(timePoints > 1)
        self.timeIndexSpinBox.setVisible(timePoints > 1)
        _matrix, source = self.logic.ijkToRasMatrix(image)
        labels = labelGroupNames(self.path)
        orientation = _("RFC-4 orientation") if source == "rfc4" else _("axes {0}").format(source.replace("-", " "))
        self.infoLabel.text = " · ".join(
            [
                info[0]["dtype"],
                _("{n} channel(s)").format(n=channels),
                _("{n} time point(s)").format(n=timePoints),
                _("labels: {names}").format(names=", ".join(labels)) if labels else _("no labels"),
                orientation,
            ]
        )
        self.levelTable.selectRow(self.logic.selectLevel(self.multiscales, self.logic.maxBytesFromSettings()))
        self.updateLevelStatus()

    @staticmethod
    def formatBytes(size):
        return f"{size / 2**30:.2f} GiB" if size >= 2**30 else f"{size / 2**20:.1f} MiB"

    def updateLevelStatus(self):
        """Bold the level the budget selects; tick the levels present in the scene."""
        if not self.path or self.multiscales is None:
            return
        recommended = self.logic.selectLevel(self.multiscales, self.logic.maxBytesFromSettings())
        loaded, regions = set(), set()
        for node in slicer.util.getNodesByClass("vtkMRMLVolumeNode"):
            if samePath(node.GetAttribute("OMEZarr.Path"), self.path) and node.GetAttribute("OMEZarr.Level"):
                target = regions if node.GetAttribute("OMEZarr.Region") else loaded
                target.add(int(node.GetAttribute("OMEZarr.Level")))
        for row in range(self.levelTable.rowCount):
            notes = []
            if row == recommended:
                notes.append(_("Selected by the memory budget"))
            if row in loaded:
                notes.append(_("Loaded in the scene"))
            if row in regions:
                notes.append(_("A region of this level is loaded"))
            levelItem = self.levelTable.item(row, 0)
            levelItem.setText(f"{row} ✓" if row in loaded or row in regions else str(row))
            for column in range(self.levelTable.columnCount):
                item = self.levelTable.item(row, column)
                font = item.font()
                font.setBold(row == recommended)
                item.setFont(font)
                if column == 0:
                    item.setToolTip(". ".join(notes))

    def onLoad(self):
        if not self.path or self.selectedLevel() < 0:
            return
        properties = {"level": self.selectedLevel()}
        if self.timeIndexSpinBox.value >= 0:
            properties["timeIndex"] = self.timeIndexSpinBox.value
        with slicer.util.tryWithErrorDisplay(_("Failed to load level"), waitCursor=True):
            slicer.util.loadNodeFromFile(self.path, "OMEZarr", properties)
        self.updateLevelStatus()

    def describeNodes(self, prefix, nodes):
        volumes = scalarVolumes(nodes)
        if not volumes:
            return prefix
        dims = volumes[0].GetImageData().GetDimensions()
        size = sum(n.GetImageData().GetActualMemorySize() for n in volumes) * 1024
        return _("{prefix}: level {level} · {x} × {y} × {z} voxels · {size}").format(
            prefix=prefix,
            level=volumes[0].GetAttribute("OMEZarr.Level"),
            x=dims[0],
            y=dims[1],
            z=dims[2],
            size=self.formatBytes(size),
        )

    def onRefine(self):
        if not self.path:
            return
        view = self.viewSelector.currentData
        try:
            with Progress(_("Refining view...")) as progress:
                nodes = self.logic.refineView(
                    self.path, view, timeIndex=max(0, self.timeIndexSpinBox.value), progress=progress
                )
            self.statusLabel.text = self.describeNodes(view, nodes)
        except InterruptedError:
            self.statusLabel.text = _("{view}: cancelled").format(view=view)
        except ValueError as e:  # expected refusals are shown in place, not in a popup
            self.statusLabel.text = f"{view}: {e}"
        except Exception as e:  # noqa: BLE001
            slicer.util.errorDisplay(_("Failed to refine view: {error}").format(error=e))
        self.updateLevelStatus()

    def onAutoRefineToggled(self, enabled):
        if not enabled:
            if self.path:
                OMEZarrLogic.stopAutoRefine(self.path)
            return
        if not self.path:
            self.autoRefineCheckBox.checked = False
            return
        with slicer.util.tryWithErrorDisplay(_("Failed to start automatic refinement")):
            OMEZarrLogic.startAutoRefine(self.path)

    def cleanup(self):
        OMEZarrLogic.stopAutoRefine()

    def onCreateRoi(self):
        """A region of interest already placed: the middle half of what the slice view shows."""
        view = self.viewSelector.currentData
        bounds = self.logic.sliceViewRasBounds(view)
        center = [(bounds[i] + bounds[i + 1]) / 2.0 for i in (0, 2, 4)]
        size = [(bounds[i + 1] - bounds[i]) / 2.0 for i in (0, 2, 4)]
        roiNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsROINode", slicer.mrmlScene.GenerateUniqueName("OME-Zarr region")
        )
        roiNode.CreateDefaultDisplayNodes()
        roiNode.SetCenter(*center)
        roiNode.SetSize(*size)
        roiNode.GetDisplayNode().SetHandlesInteractive(True)
        roiNode.GetDisplayNode().SetFillOpacity(0.1)
        self.roiSelector.setCurrentNode(roiNode)
        self.statusLabel.text = _(
            "Drag the handles of the region in the slice views to adjust it, select a level above, "
            "then click 'Load region at the selected level'."
        )
        self.updateButtons()

    def onLoadRegion(self):
        roiNode = self.roiSelector.currentNode()
        if roiNode is None or not self.path:
            return
        try:
            with Progress(_("Loading region...")) as progress:
                nodes = self.logic.loadRegion(
                    self.path,
                    roiNode,
                    level=max(0, self.selectedLevel()),
                    timeIndex=max(0, self.timeIndexSpinBox.value),
                    progress=progress,
                )
        except InterruptedError:
            self.statusLabel.text = _("Region: cancelled")
            return
        except ValueError as e:
            self.statusLabel.text = _("Region: {error}").format(error=e)
            return
        scalars = scalarVolumes(nodes)
        if scalars:
            slicer.util.setSliceViewerLayers(background=scalars[0], fit=True)
        self.statusLabel.text = self.describeNodes(roiNode.GetName(), nodes)
        self.updateLevelStatus()


#
# Tests
#


class OMEZarrTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()
        self.tempDir = slicer.util.tempDirectory("OMEZarrTest")
        OMEZarrLogic.ensureNgffZarr()
        OMEZarrLogic.clearCache()
        self.savedSettings = {
            key: qt.QSettings().value(key)
            for key in (
                Settings.MAX_BYTES,
                Settings.ORIENTATION,
                Settings.LOAD_LABELS,
                Settings.TIME_MODE,
                Settings.DISPLAY_UNITS,
                Settings.AUTO_REFINE,
                Settings.LABELS_AS_SEGMENTATION,
                Settings.STORAGE_OPTIONS,
                Settings.DETECT_LABEL_MAPS,
            )
        }
        for key in self.savedSettings:
            qt.QSettings().remove(key)

    def tearDown(self):
        import shutil

        OMEZarrLogic.stopAutoRefine()
        shutil.rmtree(self.tempDir, ignore_errors=True)
        for key, value in self.savedSettings.items():
            if value is None:
                qt.QSettings().remove(key)
            else:
                qt.QSettings().setValue(key, value)
        OMEZarrLogic.resetDisplayUnits()

    def runTest(self):
        self.setUp()
        try:
            self.test_ReaderRegistration()
            self.test_RoundTripFromSlicerVolume()
            self.test_LevelSelection()
            self.test_RegionLoading()
            self.test_MicroscopyAxes()
            self.test_Labels()
            self.test_TimeSeries()
            self.test_Writer()
            self.test_RefineView()
            self.test_AutoRefine()
            self.test_MultiViewRefine()
            self.test_LabelsAsSegmentation()
            self.test_Bioformats2raw()
            self.test_LabelMapDetection()
            self.test_RefineSkipsWhenNotFiner()
            self.test_StorageOptions()
            self.test_SegmentationWriter()
            self.test_WidgetInspectsByItself()
            self.test_ReadsKeepTheApplicationResponsive()
            self.test_Settings()
            if os.environ.get("OMEZARR_TEST_REMOTE"):
                self.test_RemoteStore()
        finally:
            self.tearDown()
        self.delayDisplay("OMEZarr tests passed")

    # helpers

    def writeMRHeadStore(self):
        import ngff_zarr
        import SampleData

        mrHead = SampleData.SampleDataLogic().downloadMRHead()
        image = OMEZarrLogic.ngffImageFromVolumeNode(mrHead, name="MRHead")
        multiscales = ngff_zarr.to_multiscales(image, scale_factors=[2, 4])
        storePath = os.path.join(self.tempDir, "MRHead.ome.zarr")
        ngff_zarr.to_ome_zarr(storePath, multiscales, overwrite=True)
        OMEZarrLogic.clearCache()
        return mrHead, storePath

    @staticmethod
    def ijkToRasArray(volumeNode):
        matrix = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(matrix)
        return slicer.util.arrayFromVTKMatrix(matrix)

    @staticmethod
    def saveAsOmeZarr(node, path):
        # slicer.util.saveNode ignores properties["fileType"] (it tests hasattr on a dict),
        # so select the writer through the IO manager as the save dialog does.
        messages = slicer.vtkMRMLMessageCollection()
        success = slicer.app.coreIOManager().saveNodes("OMEZarr", {"nodeID": node.GetID(), "fileName": path}, messages)
        if not success:
            logging.error(messages.GetAllMessagesAsString())
        return success

    # tests

    def test_ReaderRegistration(self):
        self.delayDisplay("Reader registration")
        ioManager = slicer.app.coreIOManager()
        self.assertEqual(str(ioManager.fileTypeFromDescription("OME-Zarr image")), "OMEZarr")
        self.assertIsNone(omeZarrRootFromPath(self.tempDir))
        self.assertEqual(omeZarrRootFromPath("https://example.org/a/b.ome.zarr/"), "https://example.org/a/b.ome.zarr")
        self.assertGreater(OMEZarrLogic.maxBytesFromSettings(), 0)

    def test_RoundTripFromSlicerVolume(self):
        self.delayDisplay("Round trip Slicer volume -> OME-Zarr -> Slicer")
        mrHead, storePath = self.writeMRHeadStore()
        ioManager = slicer.app.coreIOManager()
        # This is the detection path used by Add Data and drag-and-drop.
        self.assertEqual(str(ioManager.fileType(storePath)), "OMEZarr")
        self.assertEqual(str(ioManager.fileType(os.path.join(storePath, "zarr.json"))), "OMEZarr")

        loaded = slicer.util.loadNodeFromFile(storePath, "OMEZarr")
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded.GetAttribute("OMEZarr.Level"), "0")
        self.assertEqual(loaded.GetAttribute("OMEZarr.OrientationSource"), "rfc4")
        np.testing.assert_allclose(self.ijkToRasArray(loaded), self.ijkToRasArray(mrHead), atol=1e-6)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(loaded), slicer.util.arrayFromVolume(mrHead))

        # Automatic detection without an explicit file type.
        loadedAuto = slicer.util.loadNodeFromFile(storePath)
        self.assertIsNotNone(loadedAuto)
        self.assertTrue(samePath(loadedAuto.GetAttribute("OMEZarr.Path"), storePath))

    def test_LevelSelection(self):
        self.delayDisplay("Multiscale level selection by memory budget")
        mrHead, storePath = self.writeMRHeadStore()
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        self.assertEqual(len(multiscales.images), 3)
        full = OMEZarrLogic.volumeBytes(multiscales.images[0])
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, full), 0)
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, full // 4), 1)
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, full, copies=2), 1)
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, 1), 2)

        node = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"maxBytes": full // 4})
        self.assertEqual(node.GetAttribute("OMEZarr.Level"), "1")
        np.testing.assert_allclose(np.array(node.GetSpacing()), 2.0 * np.array(mrHead.GetSpacing()), rtol=1e-6)
        # A downsampled level must cover the same physical extent.
        boundsFull, boundsLevel = [0.0] * 6, [0.0] * 6
        mrHead.GetRASBounds(boundsFull)
        node.GetRASBounds(boundsLevel)
        np.testing.assert_allclose(boundsLevel, boundsFull, atol=max(mrHead.GetSpacing()) * 2.5)

    def test_RegionLoading(self):
        self.delayDisplay("Region-of-interest loading reads a sub-block")
        mrHead, storePath = self.writeMRHeadStore()
        full = slicer.util.arrayFromVolume(mrHead)
        ijkToRas = self.ijkToRasArray(mrHead)
        # Sub-block in IJK: i 40..100, j 60..120, k 20..50 (half-open).
        i0, i1, j0, j1, k0, k1 = 40, 100, 60, 120, 20, 50
        corners = np.array([[i, j, k, 1.0] for i in (i0, i1 - 1) for j in (j0, j1 - 1) for k in (k0, k1 - 1)])
        ras = (ijkToRas @ corners.T).T[:, :3]
        roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode")
        roi.SetCenter(*((ras.min(axis=0) + ras.max(axis=0)) / 2.0))
        roi.SetSize(*(ras.max(axis=0) - ras.min(axis=0)))

        nodes = OMEZarrLogic.loadRegion(storePath, roi, level=0)
        self.assertEqual(len(nodes), 1)
        region = slicer.util.arrayFromVolume(nodes[0])
        self.assertEqual(region.shape, (k1 - k0, j1 - j0, i1 - i0))
        np.testing.assert_array_equal(region, full[k0:k1, j0:j1, i0:i1])
        expectedOrigin = (ijkToRas @ np.array([i0, j0, k0, 1.0]))[:3]
        np.testing.assert_allclose(self.ijkToRasArray(nodes[0])[:3, 3], expectedOrigin, atol=1e-6)

    def writeMicroscopyStore(self, name="cells", withTime=False):
        import ngff_zarr
        from ngff_zarr import Omero, OmeroChannel, OmeroWindow

        rng = np.random.default_rng(0)
        shape = (3, 2, 12, 30, 40) if withTime else (2, 12, 30, 40)
        dims = ("t", "c", "z", "y", "x") if withTime else ("c", "z", "y", "x")
        data = rng.integers(0, 4000, size=shape, dtype=np.uint16)
        scale = {"z": 2.0, "y": 0.25, "x": 0.25}
        translation = {"z": 10.0, "y": -5.0, "x": 3.0}
        units = {"z": "micrometer", "y": "micrometer", "x": "micrometer"}
        if withTime:
            scale["t"], translation["t"], units["t"] = 0.5, 0.0, "second"
        image = ngff_zarr.to_ngff_image(data, dims=dims, scale=scale, translation=translation, axes_units=units)
        multiscales = ngff_zarr.to_multiscales(image, scale_factors=[2])
        multiscales.metadata.omero = Omero(
            channels=[
                OmeroChannel(color="0000FF", window=OmeroWindow(0, 4000, 100, 3000), label="DAPI"),
                OmeroChannel(color="00FF00", window=OmeroWindow(0, 4000, 200, 2500), label="GFP"),
            ]
        )
        storePath = os.path.join(self.tempDir, f"{name}.ome.zarr")
        ngff_zarr.to_ome_zarr(storePath, multiscales, overwrite=True)
        OMEZarrLogic.clearCache()
        return data, storePath

    def writeLabelGroup(self, storePath, name, shape):
        """Add a two-label ``labels/<name>`` group, with image-label metadata, to a microscopy store."""
        import ngff_zarr

        labelData = np.zeros(shape, dtype=np.uint16)
        labelData[2:6, 5:15, 10:20] = 1
        labelData[6:10, 15:25, 20:35] = 3
        labelImage = ngff_zarr.to_ngff_image(
            labelData,
            dims=("z", "y", "x"),
            scale={"z": 2.0, "y": 0.25, "x": 0.25},
            translation={"z": 10.0, "y": -5.0, "x": 3.0},
            axes_units={"z": "micrometer", "y": "micrometer", "x": "micrometer"},
        )
        labelPath = os.path.join(storePath, "labels", name)
        ngff_zarr.to_ome_zarr(
            labelPath,
            ngff_zarr.to_multiscales(labelImage, scale_factors=[2], method=ngff_zarr.Methods.ITKWASM_LABEL_IMAGE),
        )
        OMEZarrLogic.addImageLabelMetadata(
            labelPath,
            {
                "version": "0.5",
                "colors": [{"label-value": 1, "rgba": [255, 0, 0, 255]}, {"label-value": 3, "rgba": [0, 0, 255, 128]}],
                "properties": [{"label-value": 1, "name": "nucleus"}, {"label-value": 3, "name": "cytoplasm"}],
            },
        )
        OMEZarrLogic.registerLabelInParent(labelPath)
        OMEZarrLogic.clearCache()
        return labelData, labelPath

    def test_MicroscopyAxes(self):
        self.delayDisplay("c,z,y,x microscopy store in micrometres, two channels")
        data, storePath = self.writeMicroscopyStore()
        nodes = OMEZarrLogic.loadImage(storePath, level=0)
        self.assertEqual([n.GetName() for n in nodes], ["cells_DAPI", "cells_GFP"])
        self.assertEqual(nodes[0].GetDisplayNode().GetColorNodeID(), "vtkMRMLColorTableNodeBlue")
        self.assertEqual(nodes[1].GetDisplayNode().GetColorNodeID(), "vtkMRMLColorTableNodeGreen")
        self.assertAlmostEqual(nodes[1].GetDisplayNode().GetWindow(), 2300.0)
        self.assertEqual(nodes[0].GetAttribute("OMEZarr.OrientationSource"), "assumed-LPS")
        for channel, node in enumerate(nodes):
            np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), data[channel])
            np.testing.assert_allclose(node.GetSpacing(), [0.25e-3, 0.25e-3, 2.0e-3])
            # LPS origin (3, -5, 10) um -> RAS (-3, 5, 10) um -> mm
            np.testing.assert_allclose(node.GetOrigin(), [-3.0e-3, 5.0e-3, 10.0e-3])
            self.assertEqual(node.GetAttribute("OMEZarr.LengthUnit"), "micrometer")

        # The RAS assumption keeps x/y unflipped.
        Settings.set(Settings.ORIENTATION, "RAS")
        node = OMEZarrLogic.loadImage(storePath, level=0, channels=[0])[0]
        self.assertEqual(node.GetAttribute("OMEZarr.OrientationSource"), "assumed-RAS")
        np.testing.assert_allclose(node.GetOrigin(), [3.0e-3, -5.0e-3, 10.0e-3])
        Settings.set(Settings.ORIENTATION, "LPS")

    def test_Labels(self):
        self.delayDisplay("Labels group loads as a label map with its colour table")
        data, storePath = self.writeMicroscopyStore("labelled")
        labelData, labelPath = self.writeLabelGroup(storePath, "nuclei", data.shape[1:])
        self.assertEqual(labelGroupNames(storePath), ["nuclei"])
        self.assertTrue(isLabelStore(labelPath))

        nodes = OMEZarrLogic.loadImage(storePath, level=0)
        labelNodes = labelMaps(nodes)
        self.assertEqual(len(nodes), 3)
        self.assertEqual(len(labelNodes), 1)
        labelNode = labelNodes[0]
        self.assertEqual(labelNode.GetName(), "labelled_nuclei")
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(labelNode), labelData)
        np.testing.assert_allclose(self.ijkToRasArray(labelNode), self.ijkToRasArray(nodes[0]), atol=1e-9)
        colorNode = labelNode.GetDisplayNode().GetColorNode()
        self.assertEqual(colorNode.GetColorName(1), "nucleus")
        self.assertEqual(colorNode.GetColorName(3), "cytoplasm")
        rgba = [0.0] * 4
        colorNode.GetColor(3, rgba)
        np.testing.assert_allclose(rgba, [0.0, 0.0, 1.0, 128 / 255.0], atol=1e-6)

        # Region loading and a lower level carry the labels along.
        nodes = OMEZarrLogic.loadImage(storePath, level=1, channels=[0])
        self.assertEqual([n.GetAttribute("OMEZarr.Level") for n in nodes], ["1", "1"])
        self.assertTrue(nodes[1].IsA("vtkMRMLLabelMapVolumeNode"))
        np.testing.assert_allclose(nodes[1].GetSpacing(), nodes[0].GetSpacing())

        # A label store dropped on its own loads as a label map too.
        direct = slicer.util.loadNodeFromFile(labelPath, "OMEZarr")
        self.assertTrue(direct.IsA("vtkMRMLLabelMapVolumeNode"))

        # The setting turns labels off.
        nodes = OMEZarrLogic.loadImage(storePath, level=0, labels=False)
        self.assertEqual(len(nodes), 2)

    def test_TimeSeries(self):
        self.delayDisplay("Time axis loads as a Sequence, or as one time point")
        data, storePath = self.writeMicroscopyStore("timelapse", withTime=True)
        nodes = OMEZarrLogic.loadImage(storePath, level=0)
        self.assertEqual(len(nodes), 2)
        browsers = slicer.util.getNodesByClass("vtkMRMLSequenceBrowserNode")
        self.assertEqual(len(browsers), 1)
        browser = browsers[0]
        sequence = browser.GetMasterSequenceNode()
        self.assertEqual(sequence.GetNumberOfDataNodes(), 3)
        self.assertEqual(sequence.GetIndexUnit(), "second")
        self.assertEqual(sequence.GetNthIndexValue(2), "1")
        browser.SetSelectedItemNumber(2)
        slicer.modules.sequences.logic().UpdateProxyNodesFromSequences(browser)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(nodes[0]), data[2, 0])
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(nodes[1]), data[2, 1])

        # Refinement follows the browser's selected time point (a coarse level is shown so
        # that there is something finer to load).
        self.assertEqual(OMEZarrLogic.currentTimeIndex(storePath), 2)
        coarse = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"level": 1, "timeIndex": 0})
        slicer.util.setSliceViewerLayers(background=coarse, fit=True)
        refinedT = OMEZarrLogic.refineView(storePath, "Red")
        self.assertEqual(refinedT[0].GetAttribute("OMEZarr.TimeIndex"), "2")
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        region = OMEZarrLogic.regionFromRasBounds(multiscales.images[0], OMEZarrLogic.sliceViewRasBounds("Red"))
        expected = data[2, 0][tuple(slice(*region[d]) for d in ("z", "y", "x"))]
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(refinedT[0]), expected)

        # Budget accounts for every time point and channel.
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        single = OMEZarrLogic.volumeBytes(multiscales.images[0])
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, single * 6, copies=6), 0)
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, single * 6 - 1, copies=6), 1)

        node = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"level": 0, "timeIndex": 1})
        self.assertFalse(node.IsA("vtkMRMLSequenceNode"))
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), data[1, 0])
        self.assertEqual(node.GetAttribute("OMEZarr.TimeIndex"), "1")

        Settings.set(Settings.TIME_MODE, "index")
        node = OMEZarrLogic.loadImage(storePath, level=0, channels=[1])[0]
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), data[0, 1])
        Settings.set(Settings.TIME_MODE, "sequence")

    def test_Writer(self):
        self.delayDisplay("Scalar and label map volumes are written and read back")
        import SampleData

        mrHead = SampleData.SampleDataLogic().downloadMRHead()
        storePath = os.path.join(self.tempDir, "written.ome.zarr")
        self.assertTrue(self.saveAsOmeZarr(mrHead, storePath))
        self.assertTrue(os.path.isfile(os.path.join(storePath, "zarr.json")))
        loaded = OMEZarrLogic.loadImage(storePath, level=0)[0]
        np.testing.assert_allclose(self.ijkToRasArray(loaded), self.ijkToRasArray(mrHead), atol=1e-6)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(loaded), slicer.util.arrayFromVolume(mrHead))
        self.assertGreater(len(OMEZarrLogic.openMultiscales(storePath).images), 1)

        labelArray = np.zeros(slicer.util.arrayFromVolume(mrHead).shape, dtype=np.uint8)
        labelArray[40:80, 100:150, 100:160] = 2
        labelNode = slicer.util.addVolumeFromArray(
            labelArray, ijkToRAS=self.ijkToRasArray(mrHead), name="mask", nodeClassName="vtkMRMLLabelMapVolumeNode"
        )
        labelNode.CreateDefaultDisplayNodes()
        colorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLColorTableNode", "maskColors")
        colorNode.SetTypeToUser()
        colorNode.SetNumberOfColors(3)
        colorNode.SetColor(2, "lesion", 1.0, 0.5, 0.0, 1.0)
        labelNode.GetDisplayNode().SetAndObserveColorNodeID(colorNode.GetID())
        labelPath = os.path.join(storePath, "labels", "mask")
        self.assertTrue(self.saveAsOmeZarr(labelNode, labelPath))
        self.assertTrue(isLabelStore(labelPath))
        self.assertEqual(labelGroupNames(storePath), ["mask"])
        imageLabel = readStoreAttributes(labelPath)["image-label"]
        self.assertEqual(imageLabel["properties"], [{"label-value": 2, "name": "lesion"}])
        self.assertEqual(imageLabel["colors"][0]["rgba"], [255, 128, 0, 255])

        OMEZarrLogic.clearCache()
        nodes = OMEZarrLogic.loadImage(storePath, level=0)
        self.assertEqual(len(nodes), 2)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(nodes[1]), labelArray)
        self.assertEqual(nodes[1].GetDisplayNode().GetColorNode().GetColorName(2), "lesion")

    def test_RefineView(self):
        self.delayDisplay("Refining a slice view reloads its block at a finer level")
        mrHead, storePath = self.writeMRHeadStore()
        full = slicer.util.arrayFromVolume(mrHead)
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        budget = OMEZarrLogic.volumeBytes(multiscales.images[0]) // 4
        coarse = OMEZarrLogic.loadImage(storePath, maxBytes=budget)[0]
        self.assertEqual(coarse.GetAttribute("OMEZarr.Level"), "1")
        slicer.util.setSliceViewerLayers(background=coarse, fit=True)

        sliceNode = slicer.app.layoutManager().sliceWidget("Red").mrmlSliceNode()
        sliceNode.SetOrientationToAxial()
        center = np.array(mrHead.GetOrigin())
        bounds = [0.0] * 6
        mrHead.GetRASBounds(bounds)
        center = np.array([(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2])
        sliceNode.JumpSliceByCentering(*center)
        sliceNode.SetFieldOfView(40.0, 40.0, 1.0)

        nodes = OMEZarrLogic.refineView(storePath, "Red", maxBytes=budget)
        self.assertEqual(len(nodes), 1)
        refined = nodes[0]
        self.assertEqual(refined.GetAttribute("OMEZarr.Level"), "0")
        self.assertEqual(refined.GetAttribute("OMEZarr.Refined"), "1")
        array = slicer.util.arrayFromVolume(refined)
        self.assertLess(array.size, full.size // 8)
        # The refined block is a verbatim sub-block of the full-resolution image.
        rasToIjk = np.linalg.inv(self.ijkToRasArray(mrHead))
        start = np.rint((rasToIjk @ np.append(self.ijkToRasArray(refined)[:3, 3], 1.0))[:3]).astype(int)
        k, j, i = array.shape
        np.testing.assert_array_equal(
            array, full[start[2] : start[2] + k, start[1] : start[1] + j, start[0] : start[0] + i]
        )
        redLogic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
        self.assertEqual(redLogic.GetForegroundLayer().GetVolumeNode().GetID(), refined.GetID())
        self.assertEqual(refined.GetAttribute("OMEZarr.RefinedView"), "Red")

        # A second refinement of the same view replaces the first.
        nodes = OMEZarrLogic.refineView(storePath, "Red", maxBytes=budget)
        self.assertEqual(len(OMEZarrLogic.refinedNodes(storePath)), 1)

    @staticmethod
    def waitFor(condition, timeoutSeconds=10.0):
        import time

        deadline = time.time() + timeoutSeconds
        while time.time() < deadline:
            slicer.app.processEvents()
            if condition():
                return True
            time.sleep(0.05)
        return condition()

    def test_AutoRefine(self):
        self.delayDisplay("Automatic refinement follows the Red view")
        mrHead, storePath = self.writeMRHeadStore()
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        budget = OMEZarrLogic.volumeBytes(multiscales.images[0]) // 4
        coarse = OMEZarrLogic.loadImage(storePath, maxBytes=budget)[0]
        slicer.util.setSliceViewerLayers(background=coarse, fit=True)
        sliceNode = slicer.app.layoutManager().sliceWidget("Red").mrmlSliceNode()
        sliceNode.SetOrientationToAxial()
        bounds = [0.0] * 6
        mrHead.GetRASBounds(bounds)
        center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
        sliceNode.JumpSliceByCentering(*center)
        sliceNode.SetFieldOfView(40.0, 40.0, 1.0)

        def refined():
            return OMEZarrLogic.refinedNodes(storePath)

        refiner = OMEZarrLogic.startAutoRefine(storePath, "Red", delayMs=200, maxBytes=budget)
        # The layout may still settle while the first block loads; wait until the refiner rests.
        self.assertTrue(self.waitFor(lambda: refiner.refreshCount >= 1 and refiner.idle()))
        first = refined()
        self.assertEqual(len(first), 1)
        firstOrigin = np.array(first[0].GetOrigin())

        # A still view does not reload.
        settled = refiner.refreshCount
        self.assertFalse(self.waitFor(lambda: refiner.refreshCount > settled, timeoutSeconds=1.0))

        # Panning reloads once the view is still again, replacing the previous block.
        sliceNode.JumpSliceByCentering(center[0] + 25.0, center[1], center[2])
        self.assertTrue(self.waitFor(lambda: refiner.refreshCount > settled and refiner.idle()))
        second = refined()
        self.assertEqual(len(second), 1)
        self.assertGreater(np.abs(np.array(second[0].GetOrigin()) - firstOrigin).max(), 10.0)

        OMEZarrLogic.stopAutoRefine(storePath)
        self.assertIsNone(OMEZarrLogic.autoRefiner(storePath))
        stopped = refiner.refreshCount
        sliceNode.JumpSliceByCentering(center[0] - 25.0, center[1], center[2])
        self.assertFalse(self.waitFor(lambda: refiner.refreshCount > stopped, timeoutSeconds=1.0))

    def test_MultiViewRefine(self):
        self.delayDisplay("Each slice view keeps its own refined block with the coarse window/level")
        mrHead, storePath = self.writeMRHeadStore()
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        budget = OMEZarrLogic.volumeBytes(multiscales.images[0]) // 4
        coarse = OMEZarrLogic.loadImage(storePath, maxBytes=budget)[0]
        coarse.GetDisplayNode().SetAutoWindowLevel(False)
        coarse.GetDisplayNode().SetWindowLevel(123.0, 45.0)
        slicer.util.setSliceViewerLayers(background=coarse, fit=True)
        layoutManager = slicer.app.layoutManager()
        for viewName in ("Red", "Green"):
            layoutManager.sliceWidget(viewName).mrmlSliceNode().SetFieldOfView(40.0, 40.0, 1.0)
        red = OMEZarrLogic.refineView(storePath, "Red", maxBytes=budget)[0]
        green = OMEZarrLogic.refineView(storePath, "Green", maxBytes=budget)[0]
        self.assertEqual(len(OMEZarrLogic.refinedNodes(storePath)), 2)
        self.assertEqual(
            layoutManager.sliceWidget("Red").sliceLogic().GetForegroundLayer().GetVolumeNode().GetID(), red.GetID()
        )
        self.assertEqual(
            layoutManager.sliceWidget("Green").sliceLogic().GetForegroundLayer().GetVolumeNode().GetID(), green.GetID()
        )
        self.assertIsNone(layoutManager.sliceWidget("Yellow").sliceLogic().GetForegroundLayer().GetVolumeNode())
        for node in (red, green):
            self.assertAlmostEqual(node.GetDisplayNode().GetWindow(), 123.0)
            self.assertAlmostEqual(node.GetDisplayNode().GetLevel(), 45.0)
        # Refining Green again leaves Red's block alone.
        OMEZarrLogic.refineView(storePath, "Green", maxBytes=budget)
        self.assertEqual([n.GetID() for n in OMEZarrLogic.refinedNodes(storePath, "Red")], [red.GetID()])

    def test_LabelsAsSegmentation(self):
        self.delayDisplay("Labels load as a Segmentation when the setting is on")
        data, storePath = self.writeMicroscopyStore("segmented")
        _labelData, labelPath = self.writeLabelGroup(storePath, "cells", data.shape[1:])
        Settings.set(Settings.LABELS_AS_SEGMENTATION, True)
        labelMapsBefore = len(slicer.util.getNodesByClass("vtkMRMLLabelMapVolumeNode"))
        nodes = OMEZarrLogic.loadImage(storePath, level=0, channels=[0])
        Settings.set(Settings.LABELS_AS_SEGMENTATION, False)
        self.assertEqual(len(nodes), 2)
        segmentation = nodes[1]
        self.assertTrue(segmentation.IsA("vtkMRMLSegmentationNode"))
        self.assertEqual(len(slicer.util.getNodesByClass("vtkMRMLLabelMapVolumeNode")), labelMapsBefore)
        segments = segmentation.GetSegmentation()
        self.assertEqual(segments.GetNumberOfSegments(), 2)
        self.assertEqual(sorted(segments.GetNthSegment(i).GetName() for i in range(2)), ["cytoplasm", "nucleus"])
        self.assertTrue(samePath(segmentation.GetAttribute("OMEZarr.Path"), labelPath))
        exported = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", "check")
        slicer.modules.segmentations.logic().ExportAllSegmentsToLabelmapNode(
            segmentation, exported, slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY
        )
        self.assertEqual(set(np.unique(slicer.util.arrayFromVolume(exported)).tolist()), {0, 1, 2})
        slicer.mrmlScene.RemoveNode(exported)

    def test_Bioformats2raw(self):
        self.delayDisplay("A bioformats2raw container loads every image series")
        import ngff_zarr

        root = os.path.join(self.tempDir, "converted.zarr")
        os.makedirs(root)
        with open(os.path.join(root, ".zgroup"), "w", encoding="utf-8") as fp:
            json.dump({"zarr_format": 2}, fp)
        with open(os.path.join(root, ".zattrs"), "w", encoding="utf-8") as fp:
            json.dump({"bioformats2raw.layout": 3}, fp)
        arrays = []
        for index in range(2):
            rng = np.random.default_rng(index)
            array = rng.integers(0, 255, size=(4, 8 + index, 12), dtype=np.uint8)
            arrays.append(array)
            image = ngff_zarr.to_ngff_image(array, dims=("z", "y", "x"), scale={"z": 1.0, "y": 1.0, "x": 1.0})
            ngff_zarr.to_ome_zarr(
                os.path.join(root, str(index)), ngff_zarr.to_multiscales(image, scale_factors=[2]), version="0.4"
            )
        OMEZarrLogic.clearCache()

        self.assertEqual(os.path.normpath(omeZarrRootFromPath(root)), os.path.normpath(root))
        self.assertTrue(isBioformats2rawRoot(root))
        self.assertEqual(len(bioformats2rawSeries(root)), 2)
        self.assertEqual(str(slicer.app.coreIOManager().fileType(root)), "OMEZarr")
        nodes = OMEZarrLogic.loadImage(root, level=0)
        self.assertEqual([n.GetName() for n in nodes], ["converted_0", "converted_1"])
        for node, array in zip(nodes, arrays, strict=True):
            np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), array)
        first = slicer.util.loadNodeFromFile(root, "OMEZarr", {"level": 0})
        self.assertIsNotNone(first)

    def test_LabelMapDetection(self):
        self.delayDisplay("A plain integer mask store loads as a label map")
        import ngff_zarr

        mask = np.zeros((12, 30, 40), dtype=np.uint8)
        mask[3:9, 8:20, 10:30] = 1
        image = ngff_zarr.to_ngff_image(mask, dims=("z", "y", "x"), scale={"z": 2.0, "y": 0.25, "x": 0.25})
        storePath = os.path.join(self.tempDir, "Mask.ome.zarr")
        ngff_zarr.to_ome_zarr(storePath, ngff_zarr.to_multiscales(image, scale_factors=[2]))
        OMEZarrLogic.clearCache()
        node = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"level": 0})
        self.assertTrue(node.IsA("vtkMRMLLabelMapVolumeNode"))
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), mask)
        budgeted = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"maxBytes": mask.nbytes // 2})
        self.assertTrue(budgeted.IsA("vtkMRMLLabelMapVolumeNode"))
        self.assertEqual(budgeted.GetAttribute("OMEZarr.Level"), "1")
        # Explicit override, and the setting.
        node = slicer.util.loadNodeFromFile(storePath, "OMEZarr", {"level": 0, "asLabelMap": False})
        self.assertFalse(node.IsA("vtkMRMLLabelMapVolumeNode"))
        Settings.set(Settings.DETECT_LABEL_MAPS, False)
        node = OMEZarrLogic.loadImage(storePath, level=0)[0]
        self.assertFalse(node.IsA("vtkMRMLLabelMapVolumeNode"))
        Settings.set(Settings.DETECT_LABEL_MAPS, True)
        # A microscopy intensity image is not mistaken for labels.
        data, cellsPath = self.writeMicroscopyStore("intensity")
        self.assertFalse(OMEZarrLogic.looksLikeLabelMap(OMEZarrLogic.openMultiscales(cellsPath)))

    def test_RefineSkipsWhenNotFiner(self):
        self.delayDisplay("Refinement is refused when no finer level fits the budget")
        mrHead, storePath = self.writeMRHeadStore()
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        budget = OMEZarrLogic.volumeBytes(multiscales.images[0]) // 4
        coarse = OMEZarrLogic.loadImage(storePath, maxBytes=budget)[0]
        slicer.util.setSliceViewerLayers(background=coarse, fit=True)
        before = len(OMEZarrLogic.refinedNodes(storePath))
        with self.assertRaises(ValueError):
            OMEZarrLogic.refineView(storePath, "Red", maxBytes=budget // 64)
        self.assertEqual(len(OMEZarrLogic.refinedNodes(storePath)), before)

    def test_StorageOptions(self):
        self.delayDisplay("Remote storage options: configured JSON, anonymous S3 without credentials")
        self.assertIsNone(OMEZarrLogic.storageOptions("/local/store.ome.zarr"))
        saved = {name: os.environ.pop(name, None) for name in ("AWS_ACCESS_KEY_ID", "AWS_PROFILE", "AWS_SESSION_TOKEN")}
        try:
            options = OMEZarrLogic.storageOptions("s3://bucket/store.ome.zarr")
            if not os.path.isfile(os.path.expanduser("~/.aws/credentials")):
                self.assertEqual(options, {"anon": True})
            self.assertIsNone(OMEZarrLogic.storageOptions("https://host/store.ome.zarr"))
            Settings.set(Settings.STORAGE_OPTIONS, '{"region": "us-west-2", "anon": false}')
            options = OMEZarrLogic.storageOptions("s3://bucket/store.ome.zarr")
            self.assertEqual(options["region"], "us-west-2")
            self.assertFalse(options["anon"])
            # The configured options apply to every remote store.
            self.assertEqual(OMEZarrLogic.storageOptions("https://host/store.ome.zarr")["region"], "us-west-2")
        finally:
            Settings.set(Settings.STORAGE_OPTIONS, "")
            for name, value in saved.items():
                if value is not None:
                    os.environ[name] = value

    def test_SegmentationWriter(self):
        self.delayDisplay("A segmentation is written as a label store with segment names and colours")
        import SampleData

        mrHead = SampleData.SampleDataLogic().downloadMRHead()
        labelArray = np.zeros(slicer.util.arrayFromVolume(mrHead).shape, dtype=np.uint8)
        labelArray[40:80, 100:150, 100:160] = 1
        labelArray[20:30, 50:70, 60:90] = 2
        labelNode = slicer.util.addVolumeFromArray(
            labelArray, ijkToRAS=self.ijkToRasArray(mrHead), name="seed", nodeClassName="vtkMRMLLabelMapVolumeNode"
        )
        segmentation = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "brain")
        segmentation.SetReferenceImageGeometryParameterFromVolumeNode(mrHead)
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(labelNode, segmentation)
        slicer.mrmlScene.RemoveNode(labelNode)
        segments = segmentation.GetSegmentation()
        segments.GetNthSegment(0).SetName("cortex")
        segments.GetNthSegment(0).SetColor(0.9, 0.1, 0.1)
        segments.GetNthSegment(1).SetName("ventricle")
        segments.GetNthSegment(1).SetColor(0.1, 0.1, 0.9)

        storePath = os.path.join(self.tempDir, "brain.ome.zarr")
        labelMapsBefore = len(slicer.util.getNodesByClass("vtkMRMLLabelMapVolumeNode"))
        self.assertTrue(self.saveAsOmeZarr(segmentation, storePath))
        self.assertTrue(isLabelStore(storePath))
        # The temporary export is removed again.
        self.assertEqual(len(slicer.util.getNodesByClass("vtkMRMLLabelMapVolumeNode")), labelMapsBefore)
        imageLabel = readStoreAttributes(storePath)["image-label"]
        self.assertEqual([p["name"] for p in imageLabel["properties"]], ["cortex", "ventricle"])
        self.assertEqual(imageLabel["colors"][1]["rgba"][:3], [26, 26, 230])

        OMEZarrLogic.clearCache()
        loaded = slicer.util.loadNodeFromFile(storePath, "OMEZarr")
        self.assertTrue(loaded.IsA("vtkMRMLLabelMapVolumeNode"))
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(loaded), labelArray)
        self.assertEqual(loaded.GetDisplayNode().GetColorNode().GetColorName(2), "ventricle")

    def test_WidgetInspectsByItself(self):
        self.delayDisplay("The module panel lists the levels without clicking Inspect")
        if slicer.util.mainWindow() is None:
            return
        mrHead, storePath = self.writeMRHeadStore()
        slicer.util.selectModule("OMEZarr")
        widget = slicer.modules.OMEZarrWidget
        widget.pathEdit.currentPath = self.tempDir  # not a store: quiet, empty
        self.assertEqual(widget.levelTable.rowCount, 0)
        self.assertFalse(widget.loadButton.enabled)
        widget.pathEdit.currentPath = storePath
        self.assertEqual(widget.levelTable.rowCount, 3)
        self.assertEqual(widget.timeIndexSpinBox.minimum, 0)
        self.assertFalse(widget.timeIndexSpinBox.isVisibleTo(widget.parent))
        self.assertTrue(widget.loadButton.enabled)
        self.assertEqual(widget.selectedLevel(), 0)  # the default budget fits the full resolution
        self.assertTrue(widget.levelTable.item(0, 1).font().bold())
        self.assertFalse(widget.levelTable.item(1, 1).font().bold())
        self.assertEqual(widget.levelTable.item(0, 1).text(), "256 × 256 × 130")
        # Entering the module with an OME-Zarr volume displayed fills the path.
        widget.pathEdit.currentPath = self.tempDir
        node = OMEZarrLogic.loadImage(storePath, level=2)[0]
        slicer.util.setSliceViewerLayers(background=node)
        widget.enter()
        self.assertTrue(samePath(widget.pathEdit.currentPath, storePath))
        self.assertEqual(widget.levelTable.rowCount, 3)
        self.assertEqual(widget.levelTable.item(2, 0).text(), "2 ✓")
        # A region of interest is created already placed, then loaded.
        slicer.app.layoutManager().sliceWidget("Red").mrmlSliceNode().SetFieldOfView(60.0, 60.0, 1.0)
        widget.onCreateRoi()
        roiNode = widget.roiSelector.currentNode()
        self.assertIsNotNone(roiNode)
        self.assertTrue(widget.loadRegionButton.enabled)
        self.assertGreater(min(roiNode.GetSize()), 0.0)

        def regionNodes():
            nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
            return [n for n in nodes if n.GetAttribute("OMEZarr.Region") and not n.GetAttribute("OMEZarr.Refined")]

        widget.levelTable.selectRow(0)
        regionsBefore = len(regionNodes())
        widget.onLoadRegion()
        self.assertEqual(len(regionNodes()), regionsBefore + 1)
        # A refusal is reported in the panel, not in a popup.
        widget.viewSelector.setCurrentIndex(widget.viewSelector.findData("Red"))
        self.assertEqual(widget.viewSelector.currentText, "Red (Axial)")
        Settings.set(Settings.MAX_BYTES, 1024)
        widget.onRefine()
        Settings.set(Settings.MAX_BYTES, 0)
        self.assertIn("Red:", widget.statusLabel.text)

    def test_ReadsKeepTheApplicationResponsive(self):
        self.delayDisplay("Qt events are processed while a volume is read")
        mrHead, storePath = self.writeMRHeadStore()
        ticks = []
        timer = qt.QTimer()
        timer.setInterval(1)
        timer.timeout.connect(lambda: ticks.append(1))
        timer.start()
        node = OMEZarrLogic.loadImage(storePath, level=0)[0]
        timer.stop()
        self.assertGreater(len(ticks), 0)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), slicer.util.arrayFromVolume(mrHead))

    def test_Settings(self):
        self.delayDisplay("Display units follow the store when enabled")
        data, storePath = self.writeMicroscopyStore("units")
        unitNode = slicer.mrmlScene.GetNodeByID("vtkMRMLUnitNodeApplicationLength")
        OMEZarrLogic.loadImage(storePath, level=0, channels=[0])
        self.assertEqual(unitNode.GetSuffix(), "mm")
        Settings.set(Settings.DISPLAY_UNITS, True)
        OMEZarrLogic.loadImage(storePath, level=0, channels=[0])
        self.assertEqual(unitNode.GetSuffix(), "µm")
        self.assertAlmostEqual(unitNode.GetDisplayCoefficient(), 1000.0)
        OMEZarrLogic.resetDisplayUnits()
        self.assertEqual(unitNode.GetSuffix(), "mm")
        Settings.set(Settings.DISPLAY_UNITS, False)

        Settings.set(Settings.MAX_BYTES, 0)
        self.assertGreater(OMEZarrLogic.maxBytesFromSettings(), FALLBACK_MAX_BYTES // 16)
        Settings.set(Settings.MAX_BYTES, 12345)
        self.assertEqual(OMEZarrLogic.maxBytesFromSettings(), 12345)

    def test_RemoteStore(self):
        self.delayDisplay("Remote HTTPS store (IDR)")
        url = "https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0062A/6001240.zarr"
        node = slicer.util.loadNodeFromFile(url, "OMEZarr", {"level": 2})
        self.assertIsNotNone(node)
        self.assertEqual(node.GetAttribute("OMEZarr.Level"), "2")
        self.assertEqual(slicer.util.arrayFromVolume(node).shape, (236, 68, 67))
