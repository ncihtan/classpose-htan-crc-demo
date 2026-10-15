# Earlier single-slide run (CoNIC only)

Before the 41-slide dual-model run, we ran ClassPose with the CoNIC model on one slide, HTA7_954 (`dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A`), in January 2026. This folder keeps the scripts and smaller figures from that work. The TLS analysis is only here; its methods are in [`docs/tls_analysis_methods.md`](../../docs/tls_analysis_methods.md).

The scripts read their inputs from the current folder. The zarr, the cell centroids and the parquet tables are in the Synapse folder `single_slide_conic_run` (see `data/synapse_manifest.tsv`), and the slide is in `reference_slide`. They are not maintained.

| Script | Status |
|---|---|
| `tls_revised_analysis.py` | Current TLS analysis: literature-based staging, matching the methods doc |
| `tls_radius_analysis.py` | Compares 50, 100 and 200 um radii |
| `hierarchical_classification_figure.py` | Combines CoNIC and CoNSeP calls; explains why we ran both models |
| `morphometric_analysis.py` | Cell shape and size from contours |
| `find_tls.py`, `tls_graph.py`, `tls_heatmap.py` | Early TLS detection, superseded |
| `tls_maturity.py`, `tls_maturity_logscale.py` | Superseded; the methods doc corrects their interpretation |
| `analyze_cell_densities.py`, `cell_density_visualizations.py`, `plot_cell_types*.py` | Superseded by `analysis/` |
