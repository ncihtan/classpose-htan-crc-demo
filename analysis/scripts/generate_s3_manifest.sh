#!/usr/bin/env bash
#
# Generate a manifest of classpose outputs from S3, enriched with HTAN metadata
# Usage: ./generate_s3_manifest.sh [--profile PROFILE] [--bucket BUCKET_PATH]
#
set -euo pipefail

# Defaults
AWS_PROFILE="${AWS_PROFILE:-}"
S3_BUCKET="${CLASSPOSE_S3_BASE:-s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
METADATA_DIR="${PROJECT_DIR}/metadata"
OUTPUT_DIR="${PROJECT_DIR}/classpose_outputs"
MANIFEST_FILE="$(dirname "$PROJECT_DIR")/pipeline/s3_manifest.tsv"
TEMP_FILE=$(mktemp)
LOOKUP_FILE=$(mktemp)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --profile)
            AWS_PROFILE="$2"
            shift 2
            ;;
        --bucket)
            S3_BUCKET="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--profile PROFILE] [--bucket S3_BUCKET_PATH]"
            echo "  --profile  AWS profile to use (default: \$AWS_PROFILE, else aws default chain)"
            echo "  --bucket   S3 bucket path (default: $S3_BUCKET)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Bucket name (no prefix) for building full URIs from listed keys
S3_BUCKET_NAME="${S3_BUCKET#s3://}"
S3_BUCKET_NAME="${S3_BUCKET_NAME%%/*}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Listing S3 bucket: $S3_BUCKET"
echo "Using AWS profile: ${AWS_PROFILE:-<default>}"
echo ""

# ============================================================================
# Build lookup table from metadata files
# ============================================================================
echo "Building HTAN metadata lookup from files_drs.csv..."

DRS_FILE="${METADATA_DIR}/files_drs.csv"
if [[ ! -f "$DRS_FILE" ]]; then
    echo "Warning: Metadata file not found: $DRS_FILE"
    echo "HTAN identifiers will not be included in manifest."
    touch "$LOOKUP_FILE"
else
    # Create lookup: UUID -> htan_data_file_id,htan_biospecimen_id,htan_participant_id
    # files_drs.csv columns: drs_uri,name,atlas_name,biospecimen,assay_name,level,data_file_id,...
    # drs_uri format: drs://nci-crdc.datacommons.io/dg.4DFC/UUID
    tail -n +2 "$DRS_FILE" | while IFS=',' read -r drs_uri name atlas_name biospecimen assay_name level data_file_id rest; do
        # Extract UUID from drs_uri (after dg.4DFC/)
        uuid=$(echo "$drs_uri" | sed 's/.*dg\.4DFC\///')
        # Convert to S3 format: dg_4DFC_UUID
        sample_id="dg_4DFC_${uuid}"
        # Extract participant ID from biospecimen (e.g., HTA7_925_9 -> HTA7_925)
        participant_id=$(echo "$biospecimen" | sed 's/_[0-9]*$//')
        # Output: sample_id|data_file_id|biospecimen_id|participant_id
        echo "${sample_id}|${data_file_id}|${biospecimen}|${participant_id}"
    done > "$LOOKUP_FILE"

    LOOKUP_COUNT=$(wc -l < "$LOOKUP_FILE" | tr -d ' ')
    echo "  Loaded $LOOKUP_COUNT sample mappings"
fi

# ============================================================================
# List S3 contents
# ============================================================================
echo ""
echo "Fetching file listing from S3..."
aws s3 ls "$S3_BUCKET/" ${AWS_PROFILE:+--profile "$AWS_PROFILE"} --recursive > "$TEMP_FILE"

# Count total files
TOTAL_FILES=$(wc -l < "$TEMP_FILE" | tr -d ' ')
echo "Found $TOTAL_FILES total entries"

# ============================================================================
# Generate manifest with HTAN metadata
# ============================================================================
echo "Generating manifest..."

# Write manifest header
echo -e "sample_id\thtan_data_file_id\thtan_biospecimen_id\thtan_participant_id\tmodel\tfile_type\tfilename\ts3_key\ts3_uri\tsize_bytes\tdate_modified" > "$MANIFEST_FILE"

# Process each line and extract metadata
while IFS= read -r line; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Parse S3 ls output: date time size key
    date_part=$(echo "$line" | awk '{print $1}')
    time_part=$(echo "$line" | awk '{print $2}')
    size=$(echo "$line" | awk '{print $3}')
    s3_key=$(echo "$line" | awk '{print $4}')

    # Skip directories (size 0 and ends with /)
    if [[ "$size" == "0" ]] || [[ -z "$s3_key" ]]; then
        continue
    fi

    # Skip .zarr internal files
    if [[ "$s3_key" == *".zarr/"* ]]; then
        continue
    fi

    # Extract sample_id (dg_4DFC_UUID format)
    if [[ "$s3_key" =~ (dg_4DFC_[A-F0-9-]+) ]]; then
        sample_id="${BASH_REMATCH[1]}"
    else
        sample_id="pipeline_info"
    fi

    # Lookup HTAN identifiers
    htan_data_file_id=""
    htan_biospecimen_id=""
    htan_participant_id=""

    if [[ -s "$LOOKUP_FILE" ]] && [[ "$sample_id" != "pipeline_info" ]]; then
        lookup_result=$(grep "^${sample_id}|" "$LOOKUP_FILE" | head -1 || echo "")
        if [[ -n "$lookup_result" ]]; then
            htan_data_file_id=$(echo "$lookup_result" | cut -d'|' -f2)
            htan_biospecimen_id=$(echo "$lookup_result" | cut -d'|' -f3)
            htan_participant_id=$(echo "$lookup_result" | cut -d'|' -f4)
        fi
    fi

    # Extract model type
    if [[ "$s3_key" == *"/conic/"* ]]; then
        model="conic"
    elif [[ "$s3_key" == *"/consep/"* ]]; then
        model="consep"
    elif [[ "$s3_key" == *"/vips_converted/"* ]]; then
        model="vips_converted"
    elif [[ "$s3_key" == *"/pipeline_info/"* ]]; then
        model="pipeline_info"
    else
        model="other"
    fi

    # Extract filename
    filename=$(basename "$s3_key")

    # Determine file type
    if [[ "$filename" == *"_cell_centroids.geojson" ]]; then
        file_type="cell_centroids"
    elif [[ "$filename" == *"_cell_contours.geojson" ]]; then
        file_type="cell_contours"
    elif [[ "$filename" == *"_cell_densities.csv" ]]; then
        file_type="cell_densities"
    elif [[ "$filename" == *"_tissue_contours.geojson" ]]; then
        file_type="tissue_contours"
    elif [[ "$filename" == *".tif" ]]; then
        file_type="image_tif"
    elif [[ "$filename" == "versions.yml" ]]; then
        file_type="versions"
    elif [[ "$filename" == *"execution_report"* ]]; then
        file_type="execution_report"
    elif [[ "$filename" == *"execution_trace"* ]]; then
        file_type="execution_trace"
    elif [[ "$filename" == *"pipeline_dag"* ]]; then
        file_type="pipeline_dag"
    else
        file_type="other"
    fi

    # Build full S3 URI
    s3_uri="s3://${S3_BUCKET_NAME}/${s3_key}"

    # Write to manifest
    echo -e "${sample_id}\t${htan_data_file_id}\t${htan_biospecimen_id}\t${htan_participant_id}\t${model}\t${file_type}\t${filename}\t${s3_key}\t${s3_uri}\t${size}\t${date_part} ${time_part}" >> "$MANIFEST_FILE"

done < "$TEMP_FILE"

# Cleanup
rm -f "$TEMP_FILE" "$LOOKUP_FILE"

# Count manifest entries (excluding header)
MANIFEST_COUNT=$(($(wc -l < "$MANIFEST_FILE") - 1))

echo ""
echo "=== Manifest Summary ==="
echo "Output file: $MANIFEST_FILE"
echo "Total entries: $MANIFEST_COUNT"
echo ""

# Print summary statistics
echo "=== By Model ==="
tail -n +2 "$MANIFEST_FILE" | cut -f5 | sort | uniq -c | sort -rn

echo ""
echo "=== By File Type ==="
tail -n +2 "$MANIFEST_FILE" | cut -f6 | sort | uniq -c | sort -rn

echo ""
echo "=== Unique Samples ==="
SAMPLE_COUNT=$(tail -n +2 "$MANIFEST_FILE" | cut -f1 | grep -v "pipeline_info" | sort -u | wc -l | tr -d ' ')
echo "$SAMPLE_COUNT unique sample IDs"

echo ""
echo "=== HTAN Participants ==="
PARTICIPANT_COUNT=$(tail -n +2 "$MANIFEST_FILE" | cut -f4 | grep -v "^$" | sort -u | wc -l | tr -d ' ')
echo "$PARTICIPANT_COUNT unique HTAN participant IDs"

echo ""
echo "=== Total Size ==="
TOTAL_BYTES=$(tail -n +2 "$MANIFEST_FILE" | awk -F'\t' '{sum += $10} END {print sum}')
TOTAL_GB=$(echo "scale=2; $TOTAL_BYTES / 1024 / 1024 / 1024" | bc)
echo "${TOTAL_GB} GB"

echo ""
echo "=== Sample Preview ==="
head -5 "$MANIFEST_FILE" | cut -f1-6 | column -t -s $'\t'

echo ""
echo "Done! Manifest saved to: $MANIFEST_FILE"
