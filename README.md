# SlicerOMEZarr

Prototype [3D Slicer](https://slicer.org) extension that opens
[OME-Zarr](https://ngff.openmicroscopy.org/) (OME-NGFF) images directly.
Reading is done by [ngff-zarr](https://github.com/fideus-labs/ngff-zarr).

I like working in Slicer, and until now I had to convert every mouse-brain
OME-Zarr dataset to NIfTI just to look at it. This extension removes that step.

## Usage

* Drag an `.ome.zarr` directory onto the Slicer window and pick
  "Load OME-Zarr image".
* `File → Add Data` and select the store's `zarr.json`.
* The **OME-Zarr** module lists the resolution levels of a store, loads a
  chosen level, or loads only the region under a Markups ROI.
* From Python:

  ```python
  slicer.util.loadNodeFromFile("/data/brain.ome.zarr", "OMEZarr", {"level": 1})
  slicer.util.loadNodeFromFile("https://host/brain.ome.zarr", "OMEZarr")
  ```

The `ngff-zarr[remote]` Python package is installed into Slicer's Python on
first use.

## What works

* Multiscales: the finest level whose volume fits the memory budget is loaded
  (1 GiB by default, adjustable in the module). A message says which level
  was chosen.
* Axes `t`, `c`, `z`, `y`, `x` in any order. Channels become separate volumes
  named, coloured and windowed from the OMERO metadata. 2D images are loaded
  as single-slice volumes.
* Spacing and origin are converted from the axis units to millimetres. RFC-4
  anatomical orientation becomes the IJK→RAS direction. Without RFC-4
  metadata, `x`/`y`/`z` are treated as LPS axes, matching ngff-zarr's ITK
  conversion, so Slicer → OME-Zarr → Slicer round trips are exact.
* Local directories, `.ozx` files, `https://` and `s3://` stores.
* Region-of-interest loading reads only the chunks that intersect the region,
  at any level.

## Still to do

* `Add Data → Choose Directory to Add` lists the chunk files instead of the
  store. This needs a change in Slicer core; use drag-and-drop or the
  `zarr.json` file meanwhile.
* Time axis as a Sequence node. Only one time index is loaded today.
* Writing volumes back to OME-Zarr.
* Loading in the background with a progress bar.
* Automatic level and region selection from the current view, for datasets
  that do not fit in memory.
* ngff-zarr should return RFC-4 orientation and channel names on read. The
  extension reads the raw metadata for now.
* Tests on macOS and Windows, CI, icon, submission to the Extensions Index.

## Development

```bash
Slicer --additional-module-paths /path/to/SlicerOMEZarr/OMEZarr
Testing/run_headless_test.sh /path/to/Slicer                       # module self-test under Xvfb
OMEZARR_TEST_REMOTE=1 Testing/run_headless_test.sh /path/to/Slicer # also test an IDR HTTPS store
```
