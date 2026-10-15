"""
Submit a Vertex AI Custom Training Job for PanNuke × Classpose.

Usage:
  uv run python vertex_job.py [--epochs 500] [--dry-run]
"""

import argparse
from google.cloud.aiplatform_v1 import JobServiceClient
from google.cloud.aiplatform_v1.types import (
    CustomJob,
    CustomJobSpec,
    WorkerPoolSpec,
    MachineSpec,
    ContainerSpec,
)

PROJECT = "htan-dcc"
REGION = "us-central1"
BUCKET = "gs://htan-dcc-pannuke"
JOB_NAME = "pannuke-classpose"
TRAIN_IMAGE = "us-central1-docker.pkg.dev/htan-dcc/pannuke/trainer:latest"


def submit(epochs: int, dry_run: bool):
    training_args = [
        "--data_dir", "/gcs/htan-dcc-pannuke/data/train",
        "--train_fraction", "1.0",
        "--epochs", str(epochs),
        "--batch_size", "8",
        "--device", "0",
        "--output_dir", "/gcs/htan-dcc-pannuke/checkpoints",
        "--model_name", f"pannuke_{epochs}ep",
        "--save_every", "50",
        "--augmentation_strategy", "enhanced",
        "--num_workers", "4",
        "--oversampling_method", "custom",
    ]

    job = CustomJob(
        display_name=f"{JOB_NAME}-{epochs}ep",
        job_spec=CustomJobSpec(
            worker_pool_specs=[
                WorkerPoolSpec(
                    machine_spec=MachineSpec(
                        machine_type="a2-highgpu-1g",
                        accelerator_type="NVIDIA_TESLA_A100",
                        accelerator_count=1,
                    ),
                    replica_count=1,
                    container_spec=ContainerSpec(
                        image_uri=TRAIN_IMAGE,
                        args=training_args,
                    ),
                )
            ],
            base_output_directory={"output_uri_prefix": f"{BUCKET}/vertex-output"},
        ),
    )

    if dry_run:
        print("Dry run — job spec:")
        print(f"  Project:  {PROJECT}")
        print(f"  Region:   {REGION}")
        print(f"  Machine:  a2-highgpu-1g + A100 (on-demand)")
        print(f"  Image:    {TRAIN_IMAGE}")
        print(f"  Epochs:   {epochs}")
        print(f"  Data:     {BUCKET}/data/train")
        print(f"  Output:   {BUCKET}/checkpoints")
        return

    client = JobServiceClient(
        client_options={"api_endpoint": f"{REGION}-aiplatform.googleapis.com"}
    )
    response = client.create_custom_job(
        parent=f"projects/{PROJECT}/locations/{REGION}",
        custom_job=job,
    )
    print(f"Job submitted: {response.name}")
    print(f"Monitor at: https://console.cloud.google.com/vertex-ai/training/custom-jobs?project={PROJECT}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    submit(args.epochs, args.dry_run)
