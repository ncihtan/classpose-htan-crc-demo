"""
Evaluate a trained Classpose model on the PanNuke fold3 test set.

Computes per-class Panoptic Quality (PQ) and compares against the PUMA baseline.

Usage:
  python evaluate.py \
    --model_path checkpoints/full/pannuke_full/pannuke_full.pth \
    --test_data_dir data/test \
    --output_dir results \
    --device mps        # or 'auto', 'cpu', '0' for CUDA

Expects data/test/ to contain images.npy and labels.npy produced by:
  python convert_pannuke.py --folds fold3 --output_dir data/test
"""

import argparse
import json
import os

import numpy as np
import torch
from tqdm import tqdm

# PanNuke classes (1-indexed as stored in labels.npy)
CLASS_NAMES = {
    0: "background",
    1: "Neoplastic",
    2: "Inflammatory",
    3: "Connective",
    4: "Dead",
    5: "Epithelial",
}

# PUMA baseline PQ scores for reference (from Classpose paper)
PUMA_BASELINE_PQ = {
    "Neoplastic": None,  # PUMA not evaluated on PanNuke classes directly
    "Inflammatory": None,
    "Connective": 0.348,  # weak stroma detection — key target for improvement
    "Dead": None,
    "Epithelial": None,
}


def match_instances(pred_inst, true_inst, iou_thresh=0.5):
    """
    Match predicted instances to ground truth instances by IoU.

    Returns:
      tp  : number of true positives
      fp  : number of false positives (unmatched predictions)
      fn  : number of false negatives (unmatched ground truth)
    """
    pred_ids = np.unique(pred_inst[pred_inst > 0])
    true_ids = np.unique(true_inst[true_inst > 0])

    if len(pred_ids) == 0 and len(true_ids) == 0:
        return 0, 0, 0

    matched_true = set()
    tp = 0

    for pid in pred_ids:
        p_mask = pred_inst == pid
        # Find overlapping GT instances
        overlapping = np.unique(true_inst[p_mask])
        overlapping = overlapping[overlapping > 0]

        best_iou = 0.0
        best_tid = None
        for tid in overlapping:
            if tid in matched_true:
                continue
            t_mask = true_inst == tid
            intersection = np.logical_and(p_mask, t_mask).sum()
            union = np.logical_or(p_mask, t_mask).sum()
            iou = intersection / union if union > 0 else 0.0
            if iou > best_iou:
                best_iou = iou
                best_tid = tid

        if best_iou >= iou_thresh and best_tid is not None:
            tp += 1
            matched_true.add(best_tid)

    fp = len(pred_ids) - tp
    fn = len(true_ids) - tp
    return tp, fp, fn


def compute_sq(pred_inst, true_inst, matched_pairs):
    """Segmentation quality: mean IoU over matched pairs."""
    if not matched_pairs:
        return 0.0
    ious = []
    for pid, tid in matched_pairs:
        p_mask = pred_inst == pid
        t_mask = true_inst == tid
        intersection = np.logical_and(p_mask, t_mask).sum()
        union = np.logical_or(p_mask, t_mask).sum()
        ious.append(intersection / union if union > 0 else 0.0)
    return float(np.mean(ious))


def match_instances_with_pairs(pred_inst, true_inst, iou_thresh=0.5):
    """Like match_instances but also returns the list of (pred_id, true_id) pairs."""
    pred_ids = np.unique(pred_inst[pred_inst > 0])
    true_ids = np.unique(true_inst[true_inst > 0])

    if len(pred_ids) == 0 and len(true_ids) == 0:
        return 0, 0, 0, []

    matched_true = set()
    pairs = []

    for pid in pred_ids:
        p_mask = pred_inst == pid
        overlapping = np.unique(true_inst[p_mask])
        overlapping = overlapping[overlapping > 0]

        best_iou = 0.0
        best_tid = None
        for tid in overlapping:
            if tid in matched_true:
                continue
            t_mask = true_inst == tid
            intersection = np.logical_and(p_mask, t_mask).sum()
            union = np.logical_or(p_mask, t_mask).sum()
            iou = intersection / union if union > 0 else 0.0
            if iou > best_iou:
                best_iou = iou
                best_tid = tid

        if best_iou >= iou_thresh and best_tid is not None:
            pairs.append((pid, best_tid))
            matched_true.add(best_tid)

    tp = len(pairs)
    fp = len(pred_ids) - tp
    fn = len(true_ids) - tp
    return tp, fp, fn, pairs


def per_class_pq(pred_inst, pred_cls, true_inst, true_cls, n_classes, iou_thresh=0.5):
    """
    Compute per-class PQ for a single image.

    pred_inst, true_inst : (H, W) instance maps
    pred_cls, true_cls   : (H, W) class maps (1-indexed, 0=bg)
    """
    results = {}
    for c in range(1, n_classes):
        # Extract instances belonging to class c
        p_c = np.zeros_like(pred_inst)
        t_c = np.zeros_like(true_inst)

        pred_ids = np.unique(pred_inst[pred_cls == c])
        pred_ids = pred_ids[pred_ids > 0]
        for pid in pred_ids:
            p_c[pred_inst == pid] = pid

        true_ids = np.unique(true_inst[true_cls == c])
        true_ids = true_ids[true_ids > 0]
        for tid in true_ids:
            t_c[true_inst == tid] = tid

        tp, fp, fn, pairs = match_instances_with_pairs(p_c, t_c, iou_thresh)
        sq = compute_sq(p_c, t_c, pairs)
        dq = tp / (tp + 0.5 * fp + 0.5 * fn) if (tp + fp + fn) > 0 else 0.0
        pq = sq * dq
        results[c] = {"tp": tp, "fp": fp, "fn": fn, "sq": sq, "dq": dq, "pq": pq}

    return results


def evaluate(model_path: str, test_data_dir: str, output_dir: str, device: str, batch_size: int, max_images: int | None = None):
    from classpose.models import ClassposeModel

    # Load test data
    images = np.load(os.path.join(test_data_dir, "images.npy"), mmap_mode="r")
    labels = np.load(os.path.join(test_data_dir, "labels.npy"), mmap_mode="r")
    n = len(images)
    if max_images is not None:
        n = min(n, max_images)
    print(f"Loaded {len(images)} test images, evaluating first {n}")

    # Detect n_classes from labels
    n_classes = int(labels[:, :, :, 1].max()) + 1
    print(f"Detected n_classes={n_classes}")

    # Resolve device
    if device == "auto":
        if torch.cuda.is_available():
            dev = torch.device("cuda")
        elif torch.backends.mps.is_available():
            dev = torch.device("mps")
        else:
            dev = torch.device("cpu")
    elif device == "mps":
        dev = torch.device("mps")
    elif device == "cpu":
        dev = torch.device("cpu")
    elif device.isdigit():
        dev = torch.device(f"cuda:{int(device)}")
    else:
        raise ValueError(f"Unknown device: {device}")

    print(f"Using device: {dev}")

    model = ClassposeModel(
        gpu=False,
        pretrained_model=model_path,
        device=dev,
        nclasses=n_classes,
    )
    model.net.eval()

    # Accumulate per-class TP/FP/FN across all images
    accum = {c: {"tp": 0, "fp": 0, "fn": 0, "sq_sum": 0.0, "n_img": 0} for c in range(1, n_classes)}

    for i in tqdm(range(n), desc="Evaluating", unit="img"):
        img = images[i].astype(np.float32)          # (H, W, 3)
        lbl = labels[i]                              # (H, W, 2)
        true_inst = lbl[:, :, 0].astype(np.int32)
        true_cls = lbl[:, :, 1].astype(np.int32)

        # Run model inference — eval() returns (masks, flows, class_masks, styles)
        pred_inst, _, pred_cls_map, _ = model.eval(img)

        pred_inst = np.array(pred_inst, dtype=np.int32)
        pred_cls_map = np.array(pred_cls_map, dtype=np.int32)

        img_results = per_class_pq(pred_inst, pred_cls_map, true_inst, true_cls, n_classes)
        for c, metrics in img_results.items():
            accum[c]["tp"] += metrics["tp"]
            accum[c]["fp"] += metrics["fp"]
            accum[c]["fn"] += metrics["fn"]
            accum[c]["sq_sum"] += metrics["sq"] * metrics["tp"]
            accum[c]["n_img"] += 1

    # Compute final per-class PQ
    final = {}
    print("\n--- Per-class PQ results ---")
    print(f"{'Class':<20} {'PQ':>6}  {'DQ':>6}  {'SQ':>6}  {'TP':>6}  {'FP':>6}  {'FN':>6}  PUMA baseline")
    for c in range(1, n_classes):
        name = CLASS_NAMES.get(c, f"class_{c}")
        tp = accum[c]["tp"]
        fp = accum[c]["fp"]
        fn = accum[c]["fn"]
        sq = accum[c]["sq_sum"] / tp if tp > 0 else 0.0
        dq = tp / (tp + 0.5 * fp + 0.5 * fn) if (tp + fp + fn) > 0 else 0.0
        pq = sq * dq
        baseline = PUMA_BASELINE_PQ.get(name, None)
        baseline_str = f"{baseline:.3f}" if baseline is not None else "  N/A"
        print(f"  {name:<18} {pq:>6.3f}  {dq:>6.3f}  {sq:>6.3f}  {tp:>6}  {fp:>6}  {fn:>6}  {baseline_str}")
        final[name] = {"pq": pq, "dq": dq, "sq": sq, "tp": tp, "fp": fp, "fn": fn}

    mean_pq = np.mean([v["pq"] for v in final.values()])
    print(f"\n  {'mPQ':<18} {mean_pq:.3f}")

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "per_class_pq.json")
    with open(out_path, "w") as f:
        json.dump({"per_class": final, "mPQ": mean_pq}, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return final, mean_pq


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate Classpose on PanNuke fold3")
    parser.add_argument("--model_path", required=True, help="Path to trained .pth checkpoint")
    parser.add_argument("--test_data_dir", default="data/test", help="Dir with fold3 images.npy and labels.npy")
    parser.add_argument("--output_dir", default="results", help="Where to write per_class_pq.json")
    parser.add_argument("--device", default="auto", help="'auto', 'cpu', 'mps', or GPU id '0'")
    parser.add_argument("--batch_size", type=int, default=1, help="Inference batch size (future use)")
    parser.add_argument("--max_images", type=int, default=None, help="Limit evaluation to first N images (for quick checks)")
    args = parser.parse_args()

    evaluate(args.model_path, args.test_data_dir, args.output_dir, args.device, args.batch_size, args.max_images)
