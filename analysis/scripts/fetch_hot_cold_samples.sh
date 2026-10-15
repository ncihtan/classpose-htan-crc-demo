#!/bin/bash
# Fetch data for immune-hot and immune-cold samples

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SAMPLE_DATA_DIR:-$PROJECT_DIR/sample_data}"
S3_BASE="${CLASSPOSE_S3_BASE:-s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model}"
mkdir -p "$OUTPUT_DIR"

echo "=============================================="
echo "Fetching immune-hot and immune-cold sample data"
echo "=============================================="

# Immune-hot: HTA7_917 (dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF)
HOT_ID="dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF"
HOT_DIR="$OUTPUT_DIR/$HOT_ID"
mkdir -p "$HOT_DIR"

echo ""
echo "Fetching immune-hot sample (HTA7_917)..."

# H&E tif
echo "  Downloading H&E tif (~2GB)..."
aws s3 cp "$S3_BASE/$HOT_ID/vips_converted/$HOT_ID.tif" \
    "$HOT_DIR/$HOT_ID.tif"

# ConSeP contours
echo "  Downloading ConSeP cell contours (~1.1GB)..."
aws s3 cp "$S3_BASE/$HOT_ID/consep/${HOT_ID}_consep_cell_contours.geojson" \
    "$HOT_DIR/${HOT_ID}_consep_cell_contours.geojson"

# ConSeP centroids
echo "  Downloading ConSeP cell centroids..."
aws s3 cp "$S3_BASE/$HOT_ID/consep/${HOT_ID}_consep_cell_centroids.geojson" \
    "$HOT_DIR/${HOT_ID}_consep_cell_centroids.geojson"

# Immune-cold: HTA7_932 (dg_4DFC_98B41CE1-C47F-4C37-BB1D-F4B49E1C4619)
COLD_ID="dg_4DFC_98B41CE1-C47F-4C37-BB1D-F4B49E1C4619"
COLD_DIR="$OUTPUT_DIR/$COLD_ID"
mkdir -p "$COLD_DIR"

echo ""
echo "Fetching immune-cold sample (HTA7_932)..."

# H&E tif
echo "  Downloading H&E tif (~2.3GB)..."
aws s3 cp "$S3_BASE/$COLD_ID/vips_converted/$COLD_ID.tif" \
    "$COLD_DIR/$COLD_ID.tif"

# ConSeP contours
echo "  Downloading ConSeP cell contours (~1.5GB)..."
aws s3 cp "$S3_BASE/$COLD_ID/consep/${COLD_ID}_consep_cell_contours.geojson" \
    "$COLD_DIR/${COLD_ID}_consep_cell_contours.geojson"

# ConSeP centroids
echo "  Downloading ConSeP cell centroids..."
aws s3 cp "$S3_BASE/$COLD_ID/consep/${COLD_ID}_consep_cell_centroids.geojson" \
    "$COLD_DIR/${COLD_ID}_consep_cell_centroids.geojson"

echo ""
echo "=============================================="
echo "Download complete!"
echo "  Hot sample: $HOT_DIR"
echo "  Cold sample: $COLD_DIR"
echo "=============================================="
