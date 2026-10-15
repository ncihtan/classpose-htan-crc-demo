"""
Convert PanNuke (HuggingFace) → Classpose .npy format.

Outputs:
  data/images.npy  shape (N, 256, 256, 3)  uint8
  data/labels.npy  shape (N, 256, 256, 2)  uint16
    channel 0: instance map  (unique int per nucleus, 0=background)
    channel 1: class map     (1-indexed class per pixel, 0=background)
                              1=Neoplastic, 2=Inflammatory, 3=Connective,
                              4=Dead, 5=Epithelial

Fold layout in the saved arrays:
  fold1: indices 0      .. n1-1
  fold2: indices n1     .. n1+n2-1
  fold3: indices n1+n2  .. N-1    ← hold-out test set

Usage:
  python convert_pannuke.py [--output_dir data] [--folds fold1 fold2 fold3]
"""

import argparse
import os

import numpy as np
from datasets import load_dataset
from tqdm import tqdm


FOLDS_DEFAULT = ["fold1", "fold2", "fold3"]
H = W = 256


def make_label(instances, categories):
    """
    Build a (H, W, 2) uint16 label from per-instance binary masks.

    instances  : list of PIL images, mode='1', shape H×W
    categories : list of int, 0-indexed PanNuke class IDs (0–4)

    Returns uint16 array (H, W, 2):
      [:, :, 0]  instance map  (1-indexed, 0=bg)
      [:, :, 1]  class map     (1-indexed, 0=bg; PanNuke class + 1)
    """
    inst = np.zeros((H, W), dtype=np.uint16)
    cls = np.zeros((H, W), dtype=np.uint16)
    for i, (mask_img, c) in enumerate(zip(instances, categories), start=1):
        m = np.array(mask_img, dtype=bool)
        inst[m] = i
        cls[m] = c + 1  # shift so 0 = background
    return np.stack([inst, cls], axis=-1)  # (H, W, 2)


def convert(output_dir: str, folds: list[str]):
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading PanNuke from HuggingFace (folds: {folds}) …")
    ds = load_dataset("RationAI/PanNuke")

    all_imgs: list[np.ndarray] = []
    all_lbls: list[np.ndarray] = []
    fold_boundaries: dict[str, tuple[int, int]] = {}

    for fold in folds:
        start = len(all_imgs)
        split = ds[fold]
        print(f"  {fold}: {len(split)} images")
        for ex in tqdm(split, desc=fold, unit="img"):
            all_imgs.append(np.array(ex["image"], dtype=np.uint8))
            all_lbls.append(make_label(ex["instances"], ex["categories"]))
        fold_boundaries[fold] = (start, len(all_imgs))

    images = np.stack(all_imgs)  # (N, 256, 256, 3)
    labels = np.stack(all_lbls)  # (N, 256, 256, 2)

    img_path = os.path.join(output_dir, "images.npy")
    lbl_path = os.path.join(output_dir, "labels.npy")
    np.save(img_path, images)
    np.save(lbl_path, labels)

    # Save fold boundaries for downstream use (evaluate.py, etc.)
    import json

    boundaries_path = os.path.join(output_dir, "fold_boundaries.json")
    with open(boundaries_path, "w") as f:
        json.dump(fold_boundaries, f, indent=2)

    print(f"\nSaved {images.shape} → {img_path}")
    print(f"Saved {labels.shape} → {lbl_path}")
    print(f"Saved fold boundaries → {boundaries_path}")
    print("\nFold boundaries (start inclusive, end exclusive):")
    for fold, (s, e) in fold_boundaries.items():
        print(f"  {fold}: [{s}, {e})  ({e - s} images)")

    # Quick sanity check
    n_bg_only = int((labels[:, :, :, 0].max(axis=(1, 2)) == 0).sum())
    print(f"\nImages with no nuclei (all-background): {n_bg_only} / {len(images)}")
    cls_map = labels[:, :, :, 1]
    print("Class pixel counts (0=bg, 1=Neopl, 2=Inflam, 3=Conn, 4=Dead, 5=Epit):")
    for c in range(6):
        n = int((cls_map == c).sum())
        print(f"  class {c}: {n:,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PanNuke to .npy for Classpose")
    parser.add_argument("--output_dir", default="data", help="Directory for output .npy files")
    parser.add_argument(
        "--folds",
        nargs="+",
        default=FOLDS_DEFAULT,
        choices=FOLDS_DEFAULT,
        help="Which folds to include (default: all three)",
    )
    args = parser.parse_args()
    convert(args.output_dir, args.folds)
