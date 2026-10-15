#!/usr/bin/env bash
#
# Fetch files from S3 based on the manifest
# Usage: ./fetch_from_s3.sh [OPTIONS]
#
# Examples:
#   ./fetch_from_s3.sh --sample dg_4DFC_06CED89D-6754-40E2-A418-4F2C288119A2
#   ./fetch_from_s3.sh --participant HTA7_925
#   ./fetch_from_s3.sh --biospecimen HTA7_925_9
#   ./fetch_from_s3.sh --model conic --file-type cell_densities --all
#   ./fetch_from_s3.sh --list-samples
#
set -euo pipefail

# Defaults
AWS_PROFILE="${AWS_PROFILE:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
MANIFEST_FILE="$(dirname "$PROJECT_DIR")/pipeline/s3_manifest.tsv"
OUTPUT_DIR="${PROJECT_DIR}/classpose_outputs/data"
DRY_RUN=false
SAMPLE_FILTER=""
HTAN_DATA_FILE_FILTER=""
HTAN_BIOSPECIMEN_FILTER=""
HTAN_PARTICIPANT_FILTER=""
MODEL_FILTER=""
FILE_TYPE_FILTER=""
FETCH_ALL=false

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Fetch classpose output files from S3 based on the manifest.

FILTER OPTIONS:
    --sample SAMPLE_ID          Filter by S3 sample ID (dg_4DFC_UUID format)
    --data-file-id ID           Filter by HTAN Data File ID (e.g., HTA7_925_1014)
    --biospecimen ID            Filter by HTAN Biospecimen ID (e.g., HTA7_925_9)
    --participant ID            Filter by HTAN Participant ID (e.g., HTA7_925)
    --model MODEL               Filter by model: conic, consep, vips_converted
    --file-type TYPE            Filter by file type: cell_centroids, cell_contours,
                                cell_densities, tissue_contours, image_tif, versions

ACTION OPTIONS:
    --all                       Fetch all matching files (required for bulk downloads)
    --output-dir DIR            Output directory (default: classpose_outputs/data)
    --dry-run                   Show what would be downloaded without downloading
    --profile PROFILE           AWS profile (default: \$AWS_PROFILE, else aws default chain)

LIST OPTIONS:
    --list-samples              List all available sample IDs with HTAN mappings
    --list-file-types           List all available file types
    --list-participants         List all HTAN participant IDs

    -h, --help                  Show this help message

EXAMPLES:
    # List available samples with HTAN IDs
    $0 --list-samples

    # Download all files for a specific HTAN participant
    $0 --participant HTA7_925 --all

    # Download cell densities for all samples (small files)
    $0 --file-type cell_densities --all

    # Download conic outputs for one biospecimen
    $0 --biospecimen HTA7_925_9 --model conic --all

    # Preview what would be downloaded
    $0 --participant HTA7_925 --model conic --dry-run --all
EOF
    exit 0
}

list_samples() {
    echo "=== Available Samples with HTAN Mappings ==="
    echo ""
    printf "%-45s %-18s %-15s %-12s\n" "SAMPLE_ID" "HTAN_DATA_FILE_ID" "BIOSPECIMEN_ID" "PARTICIPANT_ID"
    printf "%-45s %-18s %-15s %-12s\n" "---------" "-----------------" "--------------" "--------------"
    tail -n +2 "$MANIFEST_FILE" | \
        awk -F'\t' '$1 != "pipeline_info" {print $1"\t"$2"\t"$3"\t"$4}' | \
        sort -u | \
        while IFS=$'\t' read -r sample data_file bio part; do
            printf "%-45s %-18s %-15s %-12s\n" "$sample" "$data_file" "$bio" "$part"
        done
    exit 0
}

list_file_types() {
    echo "=== Available File Types ==="
    tail -n +2 "$MANIFEST_FILE" | cut -f6 | sort -u
    echo ""
    echo "=== File Type Counts ==="
    tail -n +2 "$MANIFEST_FILE" | cut -f6 | sort | uniq -c | sort -rn
    exit 0
}

list_participants() {
    echo "=== HTAN Participant IDs ==="
    tail -n +2 "$MANIFEST_FILE" | cut -f4 | grep -v "^$" | sort -u
    echo ""
    echo "=== Files per Participant ==="
    tail -n +2 "$MANIFEST_FILE" | cut -f4 | grep -v "^$" | sort | uniq -c | sort -rn
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE_FILTER="$2"
            shift 2
            ;;
        --data-file-id)
            HTAN_DATA_FILE_FILTER="$2"
            shift 2
            ;;
        --biospecimen)
            HTAN_BIOSPECIMEN_FILTER="$2"
            shift 2
            ;;
        --participant)
            HTAN_PARTICIPANT_FILTER="$2"
            shift 2
            ;;
        --model)
            MODEL_FILTER="$2"
            shift 2
            ;;
        --file-type)
            FILE_TYPE_FILTER="$2"
            shift 2
            ;;
        --all)
            FETCH_ALL=true
            shift
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --list-samples)
            list_samples
            ;;
        --list-file-types)
            list_file_types
            ;;
        --list-participants)
            list_participants
            ;;
        --profile)
            AWS_PROFILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Verify manifest exists
if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "Error: Manifest file not found: $MANIFEST_FILE"
    echo "Run ./generate_s3_manifest.sh first"
    exit 1
fi

# Build awk filter expression
# Columns: 1=sample_id, 2=htan_data_file_id, 3=htan_biospecimen_id, 4=htan_participant_id,
#          5=model, 6=file_type, 7=filename, 8=s3_key, 9=s3_uri, 10=size_bytes, 11=date_modified
AWK_FILTER="NR > 1"  # Skip header

if [[ -n "$SAMPLE_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$1 == \"$SAMPLE_FILTER\""
fi

if [[ -n "$HTAN_DATA_FILE_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$2 == \"$HTAN_DATA_FILE_FILTER\""
fi

if [[ -n "$HTAN_BIOSPECIMEN_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$3 == \"$HTAN_BIOSPECIMEN_FILTER\""
fi

if [[ -n "$HTAN_PARTICIPANT_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$4 == \"$HTAN_PARTICIPANT_FILTER\""
fi

if [[ -n "$MODEL_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$5 == \"$MODEL_FILTER\""
fi

if [[ -n "$FILE_TYPE_FILTER" ]]; then
    AWK_FILTER="$AWK_FILTER && \$6 == \"$FILE_TYPE_FILTER\""
fi

# Get matching files
MATCHING_FILES=$(awk -F'\t' "$AWK_FILTER {print \$0}" "$MANIFEST_FILE")
FILE_COUNT=$(echo "$MATCHING_FILES" | grep -c . || echo 0)

if [[ "$FILE_COUNT" -eq 0 ]]; then
    echo "No files match the specified filters."
    echo "Use --list-samples, --list-participants, or --list-file-types to see available options."
    exit 1
fi

# Calculate total size
TOTAL_BYTES=$(echo "$MATCHING_FILES" | awk -F'\t' '{sum += $10} END {print sum}')
TOTAL_MB=$(echo "scale=2; $TOTAL_BYTES / 1024 / 1024" | bc)
TOTAL_GB=$(echo "scale=2; $TOTAL_BYTES / 1024 / 1024 / 1024" | bc)

echo "=== Matching Files ==="
echo "Count: $FILE_COUNT files"
echo "Total size: ${TOTAL_MB} MB (${TOTAL_GB} GB)"
echo ""

# Show preview
echo "=== Preview (first 10) ==="
printf "%-15s %-15s %-15s %-15s %10s\n" "PARTICIPANT" "BIOSPECIMEN" "MODEL" "FILE_TYPE" "SIZE_MB"
echo "$MATCHING_FILES" | head -10 | awk -F'\t' '{printf "%-15s %-15s %-15s %-15s %10.2f\n", $4, $3, $5, $6, $10/1024/1024}'
echo ""

# Check if --all flag is required
if [[ "$FILE_COUNT" -gt 1 ]] && [[ "$FETCH_ALL" != true ]] && [[ "$DRY_RUN" != true ]]; then
    echo "Multiple files match. Use --all to download all, or refine your filters."
    echo "Use --dry-run to preview without downloading."
    exit 1
fi

if [[ "$DRY_RUN" == true ]]; then
    echo "=== Dry Run - Would Download ==="
    echo "$MATCHING_FILES" | awk -F'\t' '{print $9}'
    echo ""
    echo "Total: $FILE_COUNT files, ${TOTAL_GB} GB"
    exit 0
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=== Downloading to: $OUTPUT_DIR ==="
echo ""

# Download files
DOWNLOADED=0
FAILED=0

while IFS=$'\t' read -r sample_id htan_data_file_id htan_biospecimen_id htan_participant_id model file_type filename s3_key s3_uri size_bytes date_modified; do
    # Create subdirectory structure: data/{htan_participant_id}/{htan_biospecimen_id}/{model}/
    # Fall back to sample_id if HTAN IDs are not available
    if [[ -n "$htan_participant_id" ]]; then
        DEST_DIR="${OUTPUT_DIR}/${htan_participant_id}/${htan_biospecimen_id}/${model}"
    else
        DEST_DIR="${OUTPUT_DIR}/${sample_id}/${model}"
    fi
    mkdir -p "$DEST_DIR"
    DEST_FILE="${DEST_DIR}/${filename}"

    # Skip if file already exists and has correct size
    if [[ -f "$DEST_FILE" ]]; then
        LOCAL_SIZE=$(stat -f%z "$DEST_FILE" 2>/dev/null || stat -c%s "$DEST_FILE" 2>/dev/null || echo 0)
        if [[ "$LOCAL_SIZE" -eq "$size_bytes" ]]; then
            echo "SKIP (exists): $filename"
            ((DOWNLOADED++))
            continue
        fi
    fi

    echo "Downloading: $filename ($(echo "scale=2; $size_bytes / 1024 / 1024" | bc) MB)"

    if aws s3 cp "$s3_uri" "$DEST_FILE" ${AWS_PROFILE:+--profile "$AWS_PROFILE"} 2>/dev/null; then
        ((DOWNLOADED++))
    else
        echo "FAILED: $filename"
        ((FAILED++))
    fi

done <<< "$MATCHING_FILES"

echo ""
echo "=== Download Complete ==="
echo "Downloaded: $DOWNLOADED"
echo "Failed: $FAILED"
echo "Output directory: $OUTPUT_DIR"
