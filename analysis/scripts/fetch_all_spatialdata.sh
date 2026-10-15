#!/bin/bash
# Fetch SpatialData zarrs for all 42 samples from S3
#
# Usage: ./scripts/fetch_all_spatialdata.sh
#
# Requires: aws CLI configured with appropriate credentials

set -e

S3_BASE="${CLASSPOSE_S3_BASE:-s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SAMPLE_DATA_DIR:-$PROJECT_DIR/sample_data}"
MANIFEST="$PROJECT_DIR/data/sample_manifest.txt"

# Check if manifest exists
if [ ! -f "$MANIFEST" ]; then
    echo "Error: Sample manifest not found at $MANIFEST"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count samples
TOTAL=$(wc -l < "$MANIFEST" | tr -d ' ')
CURRENT=0

echo "Fetching SpatialData zarrs for $TOTAL samples..."
echo "Source: $S3_BASE"
echo "Destination: $OUTPUT_DIR"
echo ""

# Read sample IDs and fetch each one
while IFS= read -r SAMPLE_ID; do
    CURRENT=$((CURRENT + 1))
    DEST="$OUTPUT_DIR/$SAMPLE_ID"
    ZARR_NAME="${SAMPLE_ID}_consep_spatialdata.zarr"

    # Skip if already downloaded
    if [ -d "$DEST/$ZARR_NAME" ]; then
        echo "[$CURRENT/$TOTAL] Skipping $SAMPLE_ID (already exists)"
        continue
    fi

    echo "[$CURRENT/$TOTAL] Fetching $SAMPLE_ID..."
    mkdir -p "$DEST"

    # Download the zarr directory
    aws s3 cp --recursive --quiet \
        "$S3_BASE/$SAMPLE_ID/consep/$ZARR_NAME/" \
        "$DEST/$ZARR_NAME/"

    if [ $? -eq 0 ]; then
        echo "  -> Downloaded successfully"
    else
        echo "  -> FAILED"
    fi
done < "$MANIFEST"

echo ""
echo "Done! Downloaded zarrs to $OUTPUT_DIR/"
