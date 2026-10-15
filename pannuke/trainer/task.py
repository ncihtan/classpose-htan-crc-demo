"""
Vertex AI entrypoint. All deps are pre-installed in the Docker image.
Handles auto-resume from checkpoint on spot preemption.
"""

import os
import sys


def _find_resume_checkpoint(output_dir: str, model_name: str) -> tuple[str | None, int]:
    """
    Scan output_dir for an existing checkpoint and log file.
    Returns (checkpoint_path, completed_epochs) or (None, 0).
    """
    checkpoint = os.path.join(output_dir, model_name, model_name)
    log_path = os.path.join(output_dir, model_name, "training_run.log")

    if not os.path.exists(checkpoint):
        return None, 0

    if os.path.exists(log_path):
        last_epoch = 0
        with open(log_path) as f:
            for line in f:
                if "Epoch " in line and "Validation" not in line and "Loss:" in line:
                    try:
                        ep = int(line.split("Epoch ")[1].split(",")[0])
                        last_epoch = max(last_epoch, ep)
                    except (ValueError, IndexError):
                        pass
        if last_epoch > 0:
            print(f"  Found checkpoint at epoch {last_epoch}: {checkpoint}")
            return checkpoint, last_epoch

    return checkpoint, 0


if __name__ == "__main__":
    from run_training import main, parse_args

    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Auto-resume if a checkpoint exists (handles spot preemption)
    model_name = args.model_name or "classpose"
    checkpoint, completed = _find_resume_checkpoint(args.output_dir, model_name)
    if checkpoint and completed > 0:
        remaining = args.epochs - completed
        if remaining <= 0:
            print(f"Training already complete ({completed}/{args.epochs} epochs). Exiting.")
            sys.exit(0)
        print(f"Resuming from epoch {completed}, {remaining} epochs remaining.")
        args.resume_from = checkpoint
        args.remaining_epochs = remaining
    else:
        print("No checkpoint found, starting from scratch.")

    main(args)
