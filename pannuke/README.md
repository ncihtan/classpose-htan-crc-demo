# PanNuke training demo

This is a small test of fine-tuning ClassPose on the public PanNuke dataset (Hugging Face `RationAI/PanNuke`). The local demo trained on 2% of the data for 10 epochs, and its per-class PQ was near 0 (`results/demo/per_class_pq.json`). We tried a full run on Vertex AI in the `htan-dcc` GCP project, but none finished. Notes are in `claude.md`.

| File | What it does |
|---|---|
| `convert_pannuke.py` | Converts PanNuke to the ClassPose label format |
| `run_training.py` | Training script, adapted from ClassPose `paper_experiments` |
| `evaluate.py` | Per-class PQ on fold 3 |
| `vertex_job.py`, `trainer/`, `Dockerfile` | Vertex AI job setup |

The data is at `gs://htan-dcc-pannuke/data/`. The demo checkpoint is in the Synapse folder `pannuke`.
