# SlicerOMEZarr

[3D Slicer](https://slicer.org) extension that opens and saves
[OME-Zarr](https://ngff.openmicroscopy.org/) (OME-NGFF) images.
Reading and writing go through [ngff-zarr](https://github.com/fideus-labs/ngff-zarr).

OME-Zarr, also called OME-NGFF, is an open format for large multidimensional
images, specified by the [Open Microscopy Environment](https://www.openmicroscopy.org/)
community. An image is cut into compressed chunks, usually stored at several
resolutions, and its physical description (voxel spacing, units, orientation,
channels, labels) is kept next to the pixels as JSON. Imaging archives such as
the [Image Data Resource](https://idr.openmicroscopy.org/) and collections such
as the [OME-Zarr Open SciVis Datasets](https://github.com/InsightSoftwareConsortium/OMEZarrOpenSciVisDatasets)
publish images this way. Stores range from megabytes to terabytes and can live
on a local disk, a web server or cloud storage; the same store opens in napari,
Neuroglancer, web viewers and Python.

With this extension, in Slicer:

* **No conversion step.** Drag an `.ome.zarr` directory onto Slicer, or open an
  `https://` or `s3://` address, without making an NRRD or NIfTI copy first.
* **Images larger than memory.** The finest resolution level that fits the
  memory budget loads first. Refine a slice view, with a click or automatically
  while browsing, and it reloads what it shows at a finer level, up to full
  resolution, reading only the chunks it needs.
* **Regular Slicer data.** Images become scalar volumes, labels become label
  maps or segmentations and time series become sequences, so segmentation,
  registration and volume rendering work as usual.
* **Results other tools can read.** Volumes, label maps and segmentations are
  saved as multiscale OME-Zarr, with orientation, label names and colours.

https://github.com/user-attachments/assets/da9256d2-1c8a-4734-b562-e6e754b26ca0



**[Watch the walkthrough video](Screenshots/SlicerOMEZarr-tutorial.mp4)** (1 min, with captions), recorded on
a public CT scan of a chameleon from the
[OME-Zarr Open SciVis Datasets](https://github.com/InsightSoftwareConsortium/OMEZarrOpenSciVisDatasets),
read from S3 (2.1 GiB at 0.09 mm). Step by step with screenshots: [TUTORIAL.md](TUTORIAL.md).

## Usage

* Drag an `.ome.zarr` directory onto the Slicer window and pick
  "Load OME-Zarr image".
* `File → Add Data` and select the store's `zarr.json`.
* `File → Save`, choose the "OME-Zarr image" format for a scalar or label map
  volume. Saving a label map into `<image>.ome.zarr/labels/<name>` registers
  it as a label of that image.
* The **OME-Zarr** module lists the resolution levels of a store, loads a
  chosen level, refines what a slice view shows, or loads the region under a
  Markups ROI.
* From Python:

  ```python
  slicer.util.loadNodeFromFile("/data/brain.ome.zarr", "OMEZarr", {"level": 1})
  slicer.util.loadNodeFromFile("s3://ome-zarr-scivis/v0.5/96x2/chameleon.ome.zarr", "OMEZarr")
  ```

Requires Slicer 5.12 or newer. The `ngff-zarr[remote]` Python package, version
0.46.1 or newer, is installed into Slicer's Python on first use.

## What works

* **Multiscales**: the finest level whose volumes fit the memory budget is
  loaded. The budget defaults to a quarter of the free RAM and can be fixed in
  the module settings. A message says which level was chosen.
* **Refine view**: reloads the block shown by a slice view at the
  finest level that fits the budget, reading only the chunks it needs, and
  overlays it on the coarse volume in that view with the same window/level.
  Each slice view keeps its own block, and refinement can run automatically
  in all three views each time a view stops moving. "New ROI in view" places a
  region of interest on a view, and loading it gives an ordinary volume at the
  level you pick. Time series are refined at the
  time point selected in the sequence browser.
* **Labels**: the `labels` groups of a store load as label map volumes with
  the colours and names of their `image-label` metadata, or as Segmentation
  nodes when that setting is on. A label store can also be dropped on its own,
  and a plain integer store with few distinct values (a mask written without
  `image-label` metadata) is loaded as a label map too.
* **Time series**: the `t` axis loads as a Sequence with a browser, or as a
  single time point.
* **Axes** `t`, `c`, `z`, `y`, `x` in any order. Channels become separate
  volumes named, coloured and windowed from the OMERO metadata. 2D images are
  loaded as single-slice volumes.
* **Geometry**: spacing and origin are converted from the axis units to
  millimetres. RFC-4 anatomical orientation becomes the IJK→RAS direction.
  Without RFC-4 metadata the axes are assumed LPS (the ngff-zarr and ITK
  convention) or RAS, per the module settings, and the load log says so.
* **Display units**: optionally shows lengths in the store's unit (µm, nm).
* **Writing**: scalar volumes, label maps and segmentations as multiscale
  OME-Zarr, with RFC-4 orientation from the IJK→RAS matrix and `image-label`
  names and colours from the colour table or the segments.
* **Stores**: local directories, `.ozx` files, `https://` and `s3://`, and
  bioformats2raw containers (every image series is loaded). Public S3 buckets
  are read anonymously when no AWS credentials are configured; other options
  (region, endpoint, credentials) go in the module settings as JSON. Writing
  targets a local directory by design, as in ngff-zarr; upload it afterwards.
* **Memory**: volumes are read slab by slab straight into the VTK buffer, so a
  volume is held once in memory, not twice.
* **Progress and cancel** while reading.

## Modules

* **OME-Zarr** (Informatics): inspects the resolution levels of a store, loads
  a level or a region of interest, refines slice views at full resolution, and
  holds the settings. It also registers the OME-Zarr file reader, file writer
  and drop handler.

## Still to do

* `Add Data → Choose Directory to Add` lists the chunk files instead of the
  store, and `slicer.util.saveNode` ignores a requested file type. Both need
  changes in Slicer core (tracked in issue #1); use drag-and-drop, the
  `zarr.json` file, or the save dialog meanwhile.
* Submission to the Extensions Index once the repository is public.

## Development

```bash
Slicer --additional-module-paths /path/to/SlicerOMEZarr/OMEZarr
Testing/run_headless_test.sh /path/to/Slicer                       # module self-test under Xvfb
OMEZARR_TEST_REMOTE=1 Testing/run_headless_test.sh /path/to/Slicer # also test an IDR HTTPS store
```

The same self-test runs in GitHub Actions against the latest stable Slicer on
Linux, macOS and Windows.

## License

[MIT](LICENSE).
