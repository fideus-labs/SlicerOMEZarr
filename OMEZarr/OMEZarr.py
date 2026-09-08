"""OME-Zarr (OME-NGFF) support for 3D Slicer.

Experimental prototype. It registers:

* ``OMEZarrFileReader`` - a scripted file reader, so ``slicer.util.loadNodeFromFile``,
  ``File -> Add Data`` and drag-and-drop can open an OME-Zarr store.
* ``OMEZarrFileDialog`` - a scripted drop target that recognises dropped OME-Zarr
  directories (a store is a directory, which the generic "Add Data" dialog would
  otherwise expand into thousands of chunk files).
* ``OMEZarrLogic`` - metadata inspection, multiscale level selection, NGFF -> RAS
  geometry, and region-of-interest (chunk-aligned, partial) reading.

All NGFF parsing, multiscales handling, store access (local, zip, remote) and
RFC-4 anatomical orientation come from ``ngff-zarr``; this module only maps the
result onto MRML volume nodes.
"""

import json
import logging
import os
import re

import numpy as np
import qt
import vtk

import slicer
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

DEFAULT_MAX_BYTES = 1 << 30  # 1 GiB per loaded volume (per channel/time point)
SETTINGS_MAX_BYTES_KEY = "OMEZarr/MaxBytes"

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


#
# Module
#


class OMEZarr(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("OME-Zarr")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "Informatics")]
        self.parent.dependencies = []
        self.parent.contributors = ["Valentin Boussot (Fideus Labs)", "Matt McCormick (Fideus Labs)"]
        self.parent.helpText = _(
            "Read OME-Zarr (OME-NGFF) images. Drag an .ome.zarr directory into Slicer, or use "
            "File > Add Data. Multiscale levels are selected automatically to fit a memory budget; "
            "this module lets you inspect levels and load a region of interest at full resolution."
        )
        self.parent.acknowledgementText = _("Built on ngff-zarr (https://github.com/fideus-labs/ngff-zarr).")


#
# Store discovery helpers
#


def isRemoteUrl(path):
    return bool(re.match(r"^(https?|s3|gs|gcs|az|abfs)://", str(path)))


def _multiscalesInAttributes(attrs):
    if not isinstance(attrs, dict):
        return False
    ome = attrs.get("ome", attrs)
    return isinstance(ome, dict) and "multiscales" in ome


def _hasMultiscalesMetadata(directory):
    for name in ("zarr.json", ".zattrs"):
        metadataFile = os.path.join(directory, name)
        if not os.path.isfile(metadataFile):
            continue
        try:
            with open(metadataFile, encoding="utf-8") as fp:
                attrs = json.load(fp)
        except (OSError, ValueError):
            continue
        if name == "zarr.json":
            attrs = attrs.get("attributes", {})
        if _multiscalesInAttributes(attrs):
            return True
    return False


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
    return path if _hasMultiscalesMetadata(path) else None


#
# Logic
#


class OMEZarrLogic(ScriptedLoadableModuleLogic):
    """NGFF <-> MRML mapping. Stateless apart from a small multiscales cache."""

    _multiscalesCache = {}

    @staticmethod
    def ensureNgffZarr():
        try:
            import ngff_zarr  # noqa: F401
        except ImportError:
            if slicer.util.mainWindow() and not slicer.app.testingEnabled():
                if not slicer.util.confirmOkCancelDisplay(
                    _("The 'ngff-zarr' Python package is required to read OME-Zarr images. Install it now?")
                ):
                    raise RuntimeError("ngff-zarr is not installed") from None
            with slicer.util.tryWithErrorDisplay(_("Failed to install ngff-zarr"), waitCursor=True):
                slicer.util.pip_install("ngff-zarr[remote]")
        import ngff_zarr

        return ngff_zarr

    @classmethod
    def openMultiscales(cls, path, useCache=True):
        """Open a store with ngff-zarr; arrays stay lazy (dask), nothing is read yet."""
        ngff_zarr = cls.ensureNgffZarr()
        key = str(path)
        if useCache and key in cls._multiscalesCache:
            return cls._multiscalesCache[key]
        multiscales = ngff_zarr.from_ome_zarr(key)
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

    @classmethod
    def volumeBytes(cls, image):
        """Bytes of one spatial volume (single channel, single time point) at this level."""
        size = 1
        for n in cls.spatialShape(image).values():
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
    def selectLevel(cls, multiscales, maxBytes):
        """Finest level whose single-volume size fits the budget (coarsest level otherwise)."""
        for index, image in enumerate(multiscales.images):
            if cls.volumeBytes(image) <= maxBytes:
                return index
        return len(multiscales.images) - 1

    @staticmethod
    def maxBytesFromSettings():
        value = qt.QSettings().value(SETTINGS_MAX_BYTES_KEY)
        try:
            return int(value) if value is not None else DEFAULT_MAX_BYTES
        except (TypeError, ValueError):
            return DEFAULT_MAX_BYTES

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

    @staticmethod
    def axesOrientations(image, metadata=None):
        """RFC-4 orientation per spatial axis as AnatomicalOrientation objects.

        ``NgffImage.axes_orientations`` is used when set. Otherwise the parsed
        multiscales metadata axes are consulted (ngff-zarr's reader currently
        leaves ``axes_orientations`` unset on read, see the README).
        """
        from ngff_zarr.rfc4 import AnatomicalOrientation, AnatomicalOrientationValues

        orientations = dict(image.axes_orientations or {})
        for axis in getattr(metadata, "axes", None) or []:
            name = getattr(axis, "name", None)
            if name in orientations or name not in SPATIAL_DIMS:
                continue
            orientation = getattr(axis, "orientation", None)
            if isinstance(orientation, dict):
                value = orientation.get("value")
                if orientation.get("type", "anatomical") == "anatomical" and value:
                    try:
                        orientations[name] = AnatomicalOrientation(value=AnatomicalOrientationValues(value))
                    except ValueError:
                        logging.warning(f"Unknown anatomical orientation '{value}' for axis '{name}'")
            elif orientation is not None and hasattr(orientation, "value"):
                orientations[name] = orientation
        return orientations

    @classmethod
    def ijkToRasMatrix(cls, image, userMessages=None, metadata=None):
        """4x4 IJK->RAS (mm) for an NgffImage.

        Convention (shared with ngff-zarr's ITK conversion): NGFF x/y/z are ITK/LPS
        physical axes, ``translation`` is the origin, and the RFC-4 anatomical
        orientation of each axis (when present for all spatial axes) gives the
        direction cosines. Missing z axis (2D) gets unit spacing.
        """
        from ngff_zarr.rfc4 import anatomical_orientation_to_itk_direction

        dims = list(image.dims)
        factors = cls.unitScaleToMm(image, userMessages)
        direction = np.eye(3)
        orientations = cls.axesOrientations(image, metadata)
        spatialDims = [d for d in SPATIAL_DIMS if d in dims]
        if spatialDims and all(d in orientations for d in spatialDims):
            columns = {}
            for d in spatialDims:
                column = anatomical_orientation_to_itk_direction(orientations[d].value)
                if column is None:
                    columns = None
                    break
                columns[d] = column
            if columns:
                for columnIndex, d in enumerate(SPATIAL_DIMS):
                    if d in columns:
                        direction[:, columnIndex] = columns[d]
        spacing = [image.scale.get(d, 1.0) * factors[d] if d in dims else 1.0 for d in SPATIAL_DIMS]
        origin = [image.translation.get(d, 0.0) * factors[d] if d in dims else 0.0 for d in SPATIAL_DIMS]
        ijkToLps = np.eye(4)
        ijkToLps[:3, :3] = direction @ np.diag(spacing)
        ijkToLps[:3, 3] = origin
        lpsToRas = np.diag([-1.0, -1.0, 1.0, 1.0])
        return lpsToRas @ ijkToLps

    # ---- array access ----

    @staticmethod
    def spatialArray(image, timeIndex=0, channelIndex=0, region=None, userMessages=None):
        """Read one (z, y, x) volume as numpy. Only chunks intersecting ``region`` are read.

        ``region`` maps spatial dim -> (start, stop) in this level's index space.
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
        array = np.asarray(sub.compute() if hasattr(sub, "compute") else sub)
        if "z" not in remaining:
            array = array[np.newaxis, ...]
        return OMEZarrLogic.toVtkCompatibleDtype(array)

    @staticmethod
    def toVtkCompatibleDtype(array):
        if array.dtype == np.bool_:
            return array.astype(np.uint8)
        if array.dtype == np.float16:
            return array.astype(np.float32)
        return array

    # ---- node creation ----

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

    @classmethod
    def loadImage(
        cls,
        path,
        level=None,
        timeIndex=0,
        channels=None,
        region=None,
        name=None,
        maxBytes=None,
        userMessages=None,
        multiscales=None,
    ):
        """Create one scalar volume node per requested channel. Returns the nodes.

        ``region``: spatial dim -> (start, stop) at the chosen level (partial read).
        """
        multiscales = multiscales or cls.openMultiscales(path)
        if level is None:
            level = cls.selectLevel(multiscales, maxBytes or cls.maxBytesFromSettings())
        level = int(level)
        if level < 0 or level >= len(multiscales.images):
            raise ValueError(f"Level {level} out of range (0..{len(multiscales.images) - 1})")
        image = multiscales.images[level]
        dims = list(image.dims)

        if level > 0 and region is None and userMessages:
            full = cls.volumeBytes(multiscales.images[0])
            factor = cls.levelInfo(multiscales)[level]["downsample"]
            userMessages.AddMessage(
                vtk.vtkCommand.WarningEvent,
                f"Loaded multiscale level {level} (downsampled by "
                f"{', '.join(f'{d}: {v:g}' for d, v in factor.items())}) because level 0 "
                f"would need {full / 2**30:.2f} GiB. Use the OME-Zarr module to load a "
                "region of interest at full resolution or change the memory budget.",
            )
        if "t" in dims and image.data.shape[dims.index("t")] > 1 and userMessages:
            userMessages.AddMessage(
                vtk.vtkCommand.WarningEvent,
                f"Time axis has {image.data.shape[dims.index('t')]} points; loaded index {timeIndex}.",
            )

        ijkToRas = cls.ijkToRasMatrix(image, userMessages, multiscales.metadata)
        if region:
            start = np.array([region.get(d, (0, None))[0] for d in SPATIAL_DIMS], dtype=float)
            ijkToRas = ijkToRas.copy()
            ijkToRas[:3, 3] = (ijkToRas @ np.append(start, 1.0))[:3]

        descriptions = cls.channelDescriptions(multiscales, image)
        if channels is None:
            channels = list(range(len(descriptions)))
        baseName = name or cls.defaultNodeName(path)
        units = image.axes_units or {}
        lengthUnit = next((str(units[d]) for d in SPATIAL_DIMS if units.get(d)), "")

        nodes = []
        for channelIndex in channels:
            array = cls.spatialArray(image, timeIndex, channelIndex, region, userMessages)
            label = descriptions[channelIndex]["label"] or (f"c{channelIndex}" if len(descriptions) > 1 else None)
            nodeName = f"{baseName}_{label}" if label else baseName
            if region:
                nodeName += "_ROI"
            node = slicer.util.addVolumeFromArray(
                array, ijkToRAS=ijkToRas, name=slicer.mrmlScene.GenerateUniqueName(nodeName)
            )
            node.SetAttribute("OMEZarr.Path", str(path))
            node.SetAttribute("OMEZarr.Level", str(level))
            node.SetAttribute("OMEZarr.Dims", ",".join(dims))
            node.SetAttribute("OMEZarr.TimeIndex", str(timeIndex))
            node.SetAttribute("OMEZarr.Channel", str(channelIndex))
            node.SetAttribute("OMEZarr.LengthUnit", lengthUnit)
            if region:
                node.SetAttribute("OMEZarr.Region", ";".join(f"{d}:{s}-{e}" for d, (s, e) in region.items()))
            cls.setupDisplay(node, descriptions[channelIndex], multiChannel=len(descriptions) > 1)
            nodes.append(node)
        return nodes

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

    # ---- region of interest ----

    @classmethod
    def regionFromRasBounds(cls, image, rasBounds, metadata=None):
        """Map RAS bounds [xmin,xmax,ymin,ymax,zmin,zmax] to index ranges at this level."""
        ijkToRas = cls.ijkToRasMatrix(image, metadata=metadata)
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
    def loadRegion(cls, path, roiNode, level=0, timeIndex=0, channels=None, name=None, userMessages=None):
        multiscales = cls.openMultiscales(path)
        bounds = [0.0] * 6
        roiNode.GetRASBounds(bounds)
        region = cls.regionFromRasBounds(multiscales.images[int(level)], bounds, multiscales.metadata)
        return cls.loadImage(
            path,
            level=level,
            timeIndex=timeIndex,
            channels=channels,
            region=region,
            name=name,
            userMessages=userMessages,
            multiscales=multiscales,
        )

    # ---- write side (used by tests; candidate for an ngff-zarr helper) ----

    @staticmethod
    def ngffImageFromVolumeNode(volumeNode, name=None):
        """NgffImage (z,y,x, millimeter, RFC-4 orientation) from a scalar volume node."""
        import ngff_zarr
        from ngff_zarr.rfc4 import itk_direction_to_anatomical_orientation

        array = slicer.util.arrayFromVolume(volumeNode)
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
            array,
            dims=("z", "y", "x"),
            scale={d: float(spacing[i]) for i, d in enumerate(SPATIAL_DIMS)},
            translation={d: float(origin[i]) for i, d in enumerate(SPATIAL_DIMS)},
            name=name or volumeNode.GetName(),
            axes_units={d: "millimeter" for d in SPATIAL_DIMS},
        )
        image.axes_orientations = orientations
        return image


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
            level = properties.get("level")
            maxBytes = properties.get("maxBytes")
            with slicer.util.MessageDialog(_("Loading OME-Zarr...")) if slicer.util.mainWindow() else _NoContext():
                nodes = OMEZarrLogic.loadImage(
                    root,
                    level=int(level) if level not in (None, "") else None,
                    timeIndex=int(properties.get("timeIndex", 0)),
                    name=properties.get("name"),
                    maxBytes=int(maxBytes) if maxBytes else None,
                    userMessages=self.parent.userMessages(),
                )
        except Exception as e:  # noqa: BLE001 - report everything to the user
            import traceback

            traceback.print_exc()
            self.parent.userMessages().AddMessage(vtk.vtkCommand.ErrorEvent, f"Failed to read OME-Zarr: {e}")
            return False

        if properties.get("show", True) and nodes:
            slicer.util.setSliceViewerLayers(background=nodes[0], fit=True)
        self.parent.loadedNodes = [node.GetID() for node in nodes]
        return True


class _NoContext:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


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
# Widget (minimal: inspect levels, load a level, load a region of interest)
#


class OMEZarrWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        import ctk

        self.logic = OMEZarrLogic()
        self.multiscales = None
        self.path = None

        formLayout = qt.QFormLayout()
        self.layout.addLayout(formLayout)

        self.pathEdit = ctk.ctkPathLineEdit()
        self.pathEdit.filters = ctk.ctkPathLineEdit.Dirs
        self.pathEdit.settingKey = "OMEZarr/LastPath"
        self.pathEdit.setToolTip(_("Local .ome.zarr directory, .ozx file, or https:// / s3:// URL"))
        formLayout.addRow(_("Store:"), self.pathEdit)

        self.inspectButton = qt.QPushButton(_("Inspect"))
        formLayout.addRow(self.inspectButton)

        self.levelTable = qt.QTableWidget(0, 6)
        self.levelTable.setHorizontalHeaderLabels(
            [_("Level"), _("Shape"), _("Chunks"), _("Type"), _("Size / volume"), _("Spacing")]
        )
        self.levelTable.horizontalHeader().setStretchLastSection(True)
        self.levelTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.levelTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        formLayout.addRow(self.levelTable)

        self.levelSelector = qt.QComboBox()
        formLayout.addRow(_("Level:"), self.levelSelector)

        self.maxBytesSpinBox = qt.QSpinBox()
        self.maxBytesSpinBox.setRange(1, 1 << 20)
        self.maxBytesSpinBox.setSuffix(" MiB")
        self.maxBytesSpinBox.setValue(max(1, self.logic.maxBytesFromSettings() >> 20))
        self.maxBytesSpinBox.setToolTip(_("Automatic level selection budget used when dropping a store into Slicer"))
        formLayout.addRow(_("Memory budget:"), self.maxBytesSpinBox)

        self.timeIndexSpinBox = qt.QSpinBox()
        self.timeIndexSpinBox.setRange(0, 0)
        formLayout.addRow(_("Time index:"), self.timeIndexSpinBox)

        self.loadButton = qt.QPushButton(_("Load level"))
        formLayout.addRow(self.loadButton)

        self.roiSelector = slicer.qMRMLNodeComboBox()
        self.roiSelector.nodeTypes = ["vtkMRMLMarkupsROINode"]
        self.roiSelector.addEnabled = True
        self.roiSelector.removeEnabled = True
        self.roiSelector.noneEnabled = True
        self.roiSelector.setMRMLScene(slicer.mrmlScene)
        formLayout.addRow(_("Region of interest:"), self.roiSelector)

        self.loadRegionButton = qt.QPushButton(_("Load region at selected level"))
        self.loadRegionButton.setToolTip(_("Reads only the chunks intersecting the ROI"))
        formLayout.addRow(self.loadRegionButton)

        self.layout.addStretch(1)

        self.inspectButton.connect("clicked(bool)", self.onInspect)
        self.loadButton.connect("clicked(bool)", self.onLoad)
        self.loadRegionButton.connect("clicked(bool)", self.onLoadRegion)
        self.maxBytesSpinBox.connect("valueChanged(int)", self.onMaxBytesChanged)
        self.levelSelector.connect("currentIndexChanged(int)", self.onLevelChanged)

    def onMaxBytesChanged(self, mib):
        qt.QSettings().setValue(SETTINGS_MAX_BYTES_KEY, int(mib) << 20)

    def onInspect(self):
        path = omeZarrRootFromPath(self.pathEdit.currentPath)
        if not path:
            slicer.util.errorDisplay(_("Not an OME-Zarr multiscales store."))
            return
        with slicer.util.tryWithErrorDisplay(_("Failed to open store"), waitCursor=True):
            self.multiscales = self.logic.openMultiscales(path, useCache=False)
            self.path = path
            self.pathEdit.addCurrentPathToHistory()
        if self.multiscales is None:
            return
        info = self.logic.levelInfo(self.multiscales)
        self.levelTable.setRowCount(len(info))
        self.levelSelector.clear()
        for row, level in enumerate(info):
            values = [
                str(level["level"]),
                " x ".join(f"{d}:{n}" for d, n in zip(level["dims"], level["shape"], strict=False)),
                " x ".join(str(n) for n in level["chunks"]),
                level["dtype"],
                f"{level['bytes'] / 2**20:.1f} MiB",
                ", ".join(
                    f"{d}={level['scale'][d]:g} {level['units'].get(d) or ''}".strip()
                    for d in SPATIAL_DIMS
                    if d in level["scale"]
                ),
            ]
            for column, value in enumerate(values):
                self.levelTable.setItem(row, column, qt.QTableWidgetItem(value))
            self.levelSelector.addItem(f"{level['level']}  ({level['bytes'] / 2**20:.1f} MiB)")
        self.levelTable.resizeColumnsToContents()
        dims = list(self.multiscales.images[0].dims)
        timePoints = self.multiscales.images[0].data.shape[dims.index("t")] if "t" in dims else 1
        self.timeIndexSpinBox.setRange(0, max(0, timePoints - 1))
        self.levelSelector.setCurrentIndex(self.logic.selectLevel(self.multiscales, self.logic.maxBytesFromSettings()))

    def onLevelChanged(self, index):
        self.loadButton.setEnabled(index >= 0)

    def onLoad(self):
        if not self.path:
            self.onInspect()
            if not self.path:
                return
        with slicer.util.tryWithErrorDisplay(_("Failed to load level"), waitCursor=True):
            slicer.util.loadNodeFromFile(
                self.path,
                "OMEZarr",
                {"level": self.levelSelector.currentIndex, "timeIndex": self.timeIndexSpinBox.value},
            )

    def onLoadRegion(self):
        roiNode = self.roiSelector.currentNode()
        if roiNode is None:
            slicer.util.errorDisplay(_("Select a region of interest node first."))
            return
        if not self.path:
            self.onInspect()
            if not self.path:
                return
        with slicer.util.tryWithErrorDisplay(_("Failed to load region"), waitCursor=True):
            nodes = self.logic.loadRegion(
                self.path, roiNode, level=self.levelSelector.currentIndex, timeIndex=self.timeIndexSpinBox.value
            )
            if nodes:
                slicer.util.setSliceViewerLayers(background=nodes[0], fit=True)


#
# Tests
#


class OMEZarrTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()
        self.tempDir = slicer.util.tempDirectory("OMEZarrTest")
        OMEZarrLogic.clearCache()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tempDir, ignore_errors=True)

    def runTest(self):
        self.setUp()
        try:
            self.test_ReaderRegistration()
            self.test_RoundTripFromSlicerVolume()
            self.test_LevelSelection()
            self.test_RegionLoading()
            self.test_MicroscopyAxes()
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
        return mrHead, storePath

    @staticmethod
    def ijkToRasArray(volumeNode):
        matrix = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(matrix)
        return slicer.util.arrayFromVTKMatrix(matrix)

    # tests

    def test_ReaderRegistration(self):
        self.delayDisplay("Reader registration")
        ioManager = slicer.app.coreIOManager()
        self.assertEqual(str(ioManager.fileTypeFromDescription("OME-Zarr image")), "OMEZarr")
        self.assertIsNone(omeZarrRootFromPath(self.tempDir))
        self.assertEqual(omeZarrRootFromPath("https://example.org/a/b.ome.zarr/"), "https://example.org/a/b.ome.zarr")

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
        np.testing.assert_allclose(self.ijkToRasArray(loaded), self.ijkToRasArray(mrHead), atol=1e-6)
        np.testing.assert_array_equal(slicer.util.arrayFromVolume(loaded), slicer.util.arrayFromVolume(mrHead))

        # Automatic detection without an explicit file type.
        loadedAuto = slicer.util.loadNodeFromFile(storePath)
        self.assertIsNotNone(loadedAuto)
        self.assertEqual(loadedAuto.GetAttribute("OMEZarr.Path"), storePath)

    def test_LevelSelection(self):
        self.delayDisplay("Multiscale level selection by memory budget")
        mrHead, storePath = self.writeMRHeadStore()
        multiscales = OMEZarrLogic.openMultiscales(storePath)
        self.assertEqual(len(multiscales.images), 3)
        full = OMEZarrLogic.volumeBytes(multiscales.images[0])
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, full), 0)
        self.assertEqual(OMEZarrLogic.selectLevel(multiscales, full // 4), 1)
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

    def test_MicroscopyAxes(self):
        self.delayDisplay("c,z,y,x microscopy store in micrometres, two channels")
        import ngff_zarr

        rng = np.random.default_rng(0)
        data = rng.integers(0, 4000, size=(2, 12, 30, 40), dtype=np.uint16)
        image = ngff_zarr.to_ngff_image(
            data,
            dims=("c", "z", "y", "x"),
            scale={"z": 2.0, "y": 0.25, "x": 0.25},
            translation={"z": 10.0, "y": -5.0, "x": 3.0},
            axes_units={"z": "micrometer", "y": "micrometer", "x": "micrometer"},
        )
        multiscales = ngff_zarr.to_multiscales(image, scale_factors=[2])
        from ngff_zarr import Omero, OmeroChannel, OmeroWindow

        multiscales.metadata.omero = Omero(
            channels=[
                OmeroChannel(color="0000FF", window=OmeroWindow(0, 4000, 100, 3000), label="DAPI"),
                OmeroChannel(color="00FF00", window=OmeroWindow(0, 4000, 200, 2500), label="GFP"),
            ]
        )
        storePath = os.path.join(self.tempDir, "cells.ome.zarr")
        ngff_zarr.to_ome_zarr(storePath, multiscales, overwrite=True)

        nodes = OMEZarrLogic.loadImage(storePath, level=0)
        self.assertEqual([n.GetName() for n in nodes], ["cells_DAPI", "cells_GFP"])
        self.assertEqual(nodes[0].GetDisplayNode().GetColorNodeID(), "vtkMRMLColorTableNodeBlue")
        self.assertEqual(nodes[1].GetDisplayNode().GetColorNodeID(), "vtkMRMLColorTableNodeGreen")
        self.assertAlmostEqual(nodes[1].GetDisplayNode().GetWindow(), 2300.0)
        for channel, node in enumerate(nodes):
            np.testing.assert_array_equal(slicer.util.arrayFromVolume(node), data[channel])
            np.testing.assert_allclose(node.GetSpacing(), [0.25e-3, 0.25e-3, 2.0e-3])
            # LPS origin (3, -5, 10) um -> RAS (-3, 5, 10) um -> mm
            np.testing.assert_allclose(node.GetOrigin(), [-3.0e-3, 5.0e-3, 10.0e-3])
            self.assertEqual(node.GetAttribute("OMEZarr.LengthUnit"), "micrometer")

    def test_RemoteStore(self):
        self.delayDisplay("Remote HTTPS store (IDR)")
        url = "https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0062A/6001240.zarr"
        node = slicer.util.loadNodeFromFile(url, "OMEZarr", {"level": 2})
        self.assertIsNotNone(node)
        self.assertEqual(node.GetAttribute("OMEZarr.Level"), "2")
        self.assertEqual(slicer.util.arrayFromVolume(node).shape, (236, 68, 67))
