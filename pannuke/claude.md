# PanNuke × Classpose training exploration

Goal: train a Classpose model on PanNuke to produce a pan-cancer nucleus classifier with good coverage of tumor, immune, and stroma classes.

## Task

Explore and implement:
1. Loading PanNuke from HuggingFace
2. Converting HF format → Classpose label format
3. Setting up and running training

## Classpose

- Paper: https://www.biorxiv.org/content/10.64898/2025.12.18.695211v1
- Repo: https://github.com/sohmandal/classpose
- Training example: https://github.com/sohmandal/classpose/blob/main/paper_experiments/run_training.py
- Models on HuggingFace: https://huggingface.co/classpose/classpose

### Training API

```python
from classpose.models import ClassposeModel
from classpose.train import train_class_seg

model = ClassposeModel(gpu=True, pretrained_model="cpsam", nclasses=5)
train_class_seg(net=model.net, train_data=..., train_labels=..., n_epochs=500, ...)
```

`train_data`: list of images, shape `(H, W, 3)`
`train_labels`: list of label arrays, shape `(2, H, W)` — channel 0 is instance map (integer IDs), channel 1 is class map (class ID per pixel)

## PanNuke dataset

- HuggingFace: https://huggingface.co/datasets/RationAI/PanNuke
- 7,901 images, 256×256 px, 0.25 µm/px (40×), 189,744 nuclei, 19 cancer types, 3 folds
- 5 nucleus classes:
  - 0 = Neoplastic
  - 1 = Inflammatory
  - 2 = Connective / soft tissue (stroma)
  - 3 = Dead
  - 4 = Non-neoplastic epithelial

### Format conversion

HF stores per-instance binary masks + a categories list per image. These must be combined into Classpose's `(2, H, W)` label format:
- Channel 0: instance map — each nucleus gets a unique integer ID, background = 0
- Channel 1: class map — each pixel labeled with its nucleus class ID, background = 0

Dropping the Dead class (class 3) is optional; results in a 4-class model.

## Motivation

The existing Classpose PUMA model (closest to pan-cancer) has weak stroma detection (PQ = 0.348). PanNuke covers stroma annotations across 19 cancer types with a compatible class ontology, making it a strong training source for a pan-cancer model.

## Licence

PanNuke is CC-BY-NC-SA-4.0. The ShareAlike clause propagates to any model trained on it — released models must carry the same licence.
