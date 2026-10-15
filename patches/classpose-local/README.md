# Local ClassPose helper scripts

These three scripts sat untracked in a local clone of [sohmandal/classpose](https://github.com/sohmandal/classpose) at commit `617508c07b5c2b40be56c262c2e2a453f9811b13`. They plot ClassPose outputs for the public OpenSlide test slide `CMU-1-Small-Region.svs`. They were early experiments before nf-classpose.

| Script | What it does |
|---|---|
| `visualize_cells.py` | Overlays cell contours on the slide |
| `cell_mosaic.py` | Builds a mosaic of cell crops by class |
| `spatialdata_plots.py` | Plots from the SpatialData zarr |

Paths are relative to the current folder: put the slide in `test_data/` and the ClassPose outputs in `output/`. ClassPose itself is under CC BY-NC 4.0.
