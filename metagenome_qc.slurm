#!/bin/bash
# =============================================================================
# Metagenome Reads QC — BBDuk-based pipeline (rqcfilter2-equivalent)
#
# Replicates NMDC rqcfilter2 QC steps WITHOUT the 106 GB RQCFilterData database:
#   - Adapter/artifact trimming  (BBTools built-in adapter sequences)
#   - PhiX spike-in removal      (BBTools built-in PhiX reference)
#   - Poly-G tail trimming       (Illumina two-color chemistry artifact)
#   - Quality + length filtering
#
# Host/microbe decontamination is intentionally skipped.
#
# USAGE — submit directly or pass options via --export or env vars:
#
#   # Process a folder (auto-detects R1/R2 pairs by _R1/_R2 or _1/_2 suffix):
#   sbatch metagenome_qc.slurm --input-dir /path/to/fastq_folder
#
#   # Process a single pair explicitly:
#   sbatch metagenome_qc.slurm --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz
#
#   # Override output directory:
#   sbatch metagenome_qc.slurm --input-dir /path/to/fastq_folder --output-dir /path/to/out
#
# =============================================================================

# ---- SLURM directives -------------------------------------------------------
#SBATCH --job-name=metagenome_qc
#SBATCH --output=logs/metagenome_qc_%A_%a.out
#SBATCH --error=logs/metagenome_qc_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=msismall     
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --mail-type=ALL

# ---- Load modules -----------------------------------------------------------
module load bbmap

# ---- User configuration -----------------------------------------------------
# Set BBTOOLS_DIR to the directory containing bbduk.sh.
# MSI: bbduk.sh is on PATH via "module load bbmap" above — leave this empty.
BBTOOLS_DIR=""

# Java heap size. BBDuk only needs ~200MB for the adapter+PhiX k-mer table;
# 12g gives plenty of headroom for JVM overhead without wasting allocation.
JAVA_MEM="12g"

# Minimum read length to keep after trimming (NMDC default: 51 bp)
MIN_LEN=51

# Minimum fraction of original length to keep (NMDC default: 0.33)
MIN_LEN_FRAC=0.33

# Minimum average quality after trimming (NMDC default: 3)
MIN_AVG_QUAL=3

# Maximum Ns allowed per read (NMDC default: 3)
MAX_NS=3

# Poly-G tail trimming length threshold (NMDC default: 5)
TRIM_POLYG=5

# Number of CPU threads (should match --cpus-per-task above)
THREADS=${SLURM_CPUS_PER_TASK:-16}

# =============================================================================
# --- Do not edit below this line unless you know what you're doing -----------
# =============================================================================

set -euo pipefail

# ---- Parse arguments --------------------------------------------------------
INPUT_DIR=""
R1=""
R2=""
OUTPUT_DIR=""

print_usage() {
    echo "Usage: sbatch $0 [--input-dir DIR] [--r1 FILE --r2 FILE] [--output-dir DIR]"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)  INPUT_DIR="$2";  shift 2 ;;
        --r1)         R1="$2";         shift 2 ;;
        --r2)         R2="$2";         shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help)    print_usage ;;
        *) echo "Unknown option: $1"; print_usage ;;
    esac
done

# Validate: need either folder or explicit pair
if [[ -z "$INPUT_DIR" && ( -z "$R1" || -z "$R2" ) ]]; then
    echo "ERROR: Provide either --input-dir or both --r1 and --r2."
    print_usage
fi

# ---- Resolve BBDuk path -----------------------------------------------------
if [[ -n "$BBTOOLS_DIR" ]]; then
    BBDUK="${BBTOOLS_DIR%/}/bbduk.sh"
else
    BBDUK="bbduk.sh"
fi

if ! command -v "$BBDUK" &>/dev/null && [[ ! -x "$BBDUK" ]]; then
    echo "ERROR: bbduk.sh not found. Set BBTOOLS_DIR or ensure BBTools is on PATH."
    exit 1
fi

# ---- Build list of R1/R2 pairs ----------------------------------------------
declare -a R1_FILES R2_FILES

if [[ -n "$INPUT_DIR" ]]; then
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Input directory not found: $INPUT_DIR"
        exit 1
    fi

    # Detect R1 files by common naming conventions
    while IFS= read -r -d '' r1_path; do
        r1_base=$(basename "$r1_path")
        dir=$(dirname "$r1_path")

        # Determine matching R2 filename
        if [[ "$r1_base" =~ _R1(\.|_) ]]; then
            r2_base="${r1_base/_R1/_R2}"
        elif [[ "$r1_base" =~ _1\. ]]; then
            r2_base=$(echo "$r1_base" | sed 's/_1\./_2./')
        else
            echo "WARN: Cannot determine R2 pair for $r1_base — skipping."
            continue
        fi

        r2_path="${dir}/${r2_base}"
        if [[ ! -f "$r2_path" ]]; then
            echo "WARN: R2 file not found for $r1_base (expected $r2_path) — skipping."
            continue
        fi

        R1_FILES+=("$r1_path")
        R2_FILES+=("$r2_path")
    done < <(find "$INPUT_DIR" -maxdepth 1 \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fastq" \
                                              -o -name "*_1.fastq.gz"  -o -name "*_1.fastq"  \) \
             -print0 | sort -z)

    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No R1/R2 pairs found in $INPUT_DIR."
        exit 1
    fi

    # Default output directory next to input dir
    [[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="${INPUT_DIR%/}_qc"
else
    # Explicit pair mode
    if [[ ! -f "$R1" ]]; then echo "ERROR: R1 file not found: $R1"; exit 1; fi
    if [[ ! -f "$R2" ]]; then echo "ERROR: R2 file not found: $R2"; exit 1; fi
    R1_FILES=("$R1")
    R2_FILES=("$R2")
    [[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$(dirname "$R1")/qc_output"
fi

mkdir -p "${OUTPUT_DIR}/filtered" \
          "${OUTPUT_DIR}/stats" \
          "${OUTPUT_DIR}/histograms" \
          "${OUTPUT_DIR}/logs"

echo "======================================================================"
echo "  Metagenome QC — BBDuk"
echo "  Samples to process : ${#R1_FILES[@]}"
echo "  Output directory   : $OUTPUT_DIR"
echo "  Threads            : $THREADS"
echo "  Java heap          : $JAVA_MEM"
echo "======================================================================"

# ---- Process each sample ----------------------------------------------------
for i in "${!R1_FILES[@]}"; do
    r1="${R1_FILES[$i]}"
    r2="${R2_FILES[$i]}"

    # Derive sample name (strip path, strip extension, strip _R1/_1 suffix)
    sample=$(basename "$r1")
    sample="${sample%.fastq.gz}"
    sample="${sample%.fastq}"
    sample="${sample%_R1*}"
    sample="${sample%_1}"

    echo ""
    echo "----------------------------------------------------------------------"
    echo "  Processing sample: $sample"
    echo "  R1: $r1"
    echo "  R2: $r2"
    echo "----------------------------------------------------------------------"

    out_r1="${OUTPUT_DIR}/filtered/${sample}_R1_filtered.fastq.gz"
    out_r2="${OUTPUT_DIR}/filtered/${sample}_R2_filtered.fastq.gz"
    stats1="${OUTPUT_DIR}/stats/${sample}_filterStats.txt"
    stats2="${OUTPUT_DIR}/stats/${sample}_filterStats2.txt"
    hist_aq="${OUTPUT_DIR}/histograms/${sample}_aqhist.txt"
    hist_len="${OUTPUT_DIR}/histograms/${sample}_lhist.txt"
    tmp_r1="${OUTPUT_DIR}/filtered/${sample}_tmp_R1.fastq.gz"
    tmp_r2="${OUTPUT_DIR}/filtered/${sample}_tmp_R2.fastq.gz"

    # ------------------------------------------------------------------
    # Step 1: Adapter trimming + PhiX removal + poly-G trimming
    #   - ref=adapters  : BBTools built-in Illumina adapter database
    #   - ref=phix      : BBTools built-in PhiX174 genome
    #   - ktrim=r       : trim adapters from right end of read
    #   - k=23 mink=11  : kmer sizes for adapter detection (short overhangs)
    #   - hdist=1       : allow 1 mismatch in kmer matching
    #   - trimpolyg     : trim poly-G tails (NextSeq/NovaSeq artifact)
    # ------------------------------------------------------------------
    "$BBDUK" \
        -Xmx${JAVA_MEM} \
        in1="$r1" \
        in2="$r2" \
        out1="$tmp_r1" \
        out2="$tmp_r2" \
        ref=adapters,phix \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        trimpolyg=${TRIM_POLYG} \
        tpe tbo \
        stats="$stats1" \
        threads=${THREADS} \
        overwrite=true \
        2>&1 | tee "${OUTPUT_DIR}/stats/${sample}_step1.log"

    # ------------------------------------------------------------------
    # Step 2: Quality trimming + length + N-content filtering
    #   - qtrim=r trimq=0  : trim right end at q0 (removes very low-quality
    #                        tails; NMDC uses trimq=0 per their WDL workflow)
    #   - maq=3            : discard reads with mean quality < 3
    #   - maxns=3          : discard reads with > 3 Ns
    #   - minlen=51        : discard reads shorter than 51 bp
    #   - mlf=0.33         : OR discard if < 33% of original length remains
    #   - aqhist / lhist   : quality and length distribution stats
    # ------------------------------------------------------------------
    "$BBDUK" \
        -Xmx${JAVA_MEM} \
        in1="$tmp_r1" \
        in2="$tmp_r2" \
        out1="$out_r1" \
        out2="$out_r2" \
        qtrim=r \
        trimq=0 \
        maq=${MIN_AVG_QUAL} \
        maxns=${MAX_NS} \
        minlen=${MIN_LEN} \
        mlf=${MIN_LEN_FRAC} \
        aqhist="$hist_aq" \
        lhist="$hist_len" \
        stats="$stats2" \
        threads=${THREADS} \
        overwrite=true \
        2>&1 | tee "${OUTPUT_DIR}/stats/${sample}_step2.log"

    # Clean up intermediate files
    rm -f "$tmp_r1" "$tmp_r2"

    echo ""
    echo "  Done: $sample"
    echo "  Filtered reads : ${OUTPUT_DIR}/filtered/"
    echo "  Stats          : ${OUTPUT_DIR}/stats/"
    echo "  Histograms     : ${OUTPUT_DIR}/histograms/"
done

echo ""
echo "======================================================================"
echo "  All samples complete."
echo "======================================================================"
