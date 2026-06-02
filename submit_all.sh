#!/bin/bash
# =============================================================================
# submit_all.sh — Submit one metagenome QC job per sample pair
#
# Usage:
#   bash submit_all.sh --input-dir /path/to/reads --output-dir /path/to/out
#
# Each sample gets its own SLURM job. All jobs run in parallel (subject to
# cluster availability). Monitor with: squeue -u $USER
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/metagenome_qc.slurm"

INPUT_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)  INPUT_DIR="$2";  shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash $0 --input-dir DIR --output-dir DIR"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --input-dir and --output-dir are required."
    echo "Usage: bash $0 --input-dir DIR --output-dir DIR"
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [[ ! -f "$SLURM_SCRIPT" ]]; then
    echo "ERROR: SLURM script not found: $SLURM_SCRIPT"
    exit 1
fi

# Create output subdirectories before any jobs start (avoids race conditions
# when 40 jobs all try to mkdir simultaneously)
mkdir -p "${OUTPUT_DIR}/filtered" \
         "${OUTPUT_DIR}/stats" \
         "${OUTPUT_DIR}/histograms" \
         "${OUTPUT_DIR}/logs"

echo "======================================================================"
echo "  Submitting metagenome QC jobs"
echo "  Input  : $INPUT_DIR"
echo "  Output : $OUTPUT_DIR"
echo "======================================================================"

count=0
skipped=0

for r1 in "${INPUT_DIR}"/*_R1_001.fastq.gz; do
    # Guard against empty glob
    [[ -f "$r1" ]] || { echo "ERROR: No *_R1_001.fastq.gz files found in $INPUT_DIR"; exit 1; }

    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    if [[ ! -f "$r2" ]]; then
        echo "  WARN: No R2 found for $(basename "$r1") — skipping."
        ((skipped++))
        continue
    fi

    sample=$(basename "$r1" _R1_001.fastq.gz)

    # Name the SLURM log files after the sample for easy debugging
    job_id=$(sbatch \
        --output="${OUTPUT_DIR}/logs/${sample}_%j.out" \
        --error="${OUTPUT_DIR}/logs/${sample}_%j.err" \
        "$SLURM_SCRIPT" \
        --r1 "$r1" \
        --r2 "$r2" \
        --output-dir "$OUTPUT_DIR" \
        | awk '{print $NF}')

    echo "  Submitted $sample (job $job_id)"
    ((count++))
done

echo "======================================================================"
echo "  $count jobs submitted, $skipped skipped."
echo ""
echo "  Monitor:  squeue -u $USER"
echo "  Cancel all: scancel -u $USER"
echo "  Output:   $OUTPUT_DIR"
echo "======================================================================"
