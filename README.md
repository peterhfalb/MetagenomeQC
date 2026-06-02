# MetagenomeQC

Shotgun metagenome read quality control pipeline for SLURM-based HPC clusters. Implements the same QC steps as the [NMDC Metagenome Reads QC workflow](https://docs.microbiomedata.org/workflows/chapters/3_Metagenome_Reads_QC/index.html) using BBDuk (part of BBTools), without requiring the 106 GB RQCFilterData decontamination database.

---

## Background

Shotgun metagenomics sequences total DNA extracted from environmental or host-associated samples, enabling culture-independent characterization of microbial communities — their taxonomic composition, metabolic potential, and functional gene content. Before any downstream analysis (assembly, binning, taxonomic classification, functional annotation), raw Illumina reads must be quality-controlled to remove:

- **Sequencing adapters** — synthetic oligonucleotides ligated during library preparation that appear at read ends if the insert is shorter than the read length. Leaving adapters in reads corrupts k-mer-based analyses, alignments, and assemblies.
- **PhiX spike-in** — a bacteriophage genome added to Illumina sequencing runs as an internal quality control standard. PhiX reads are genuine sequencing data but not biological signal and must be removed.
- **Poly-G tails** — an artifact of Illumina two-color chemistry (NextSeq, NovaSeq) where signal loss at the end of a read is miscalled as repeated G bases.
- **Low-quality reads** — reads with poor average base quality or excessive ambiguous bases (N) that inflate error rates in downstream analyses.
- **Short reads** — reads too short after trimming to map reliably or contribute to assembly.

The [JGI/NMDC standard pipeline](https://docs.microbiomedata.org/workflows/chapters/3_Metagenome_Reads_QC/index.html) uses `rqcfilter2`, a wrapper around BBDuk that additionally removes host (human, mouse, dog, cat) and microbial contaminant sequences using a curated 106 GB reference database. This pipeline replicates all BBDuk QC steps but omits host/microbe decontamination, which is appropriate for environmental metagenomes where the reference database adds little value and requires substantial infrastructure.

---

## Pipeline overview

```
Raw paired-end FASTQ (R1 + R2)
         │
         ▼
┌─────────────────────────────────────────────────┐
│  Step 1 — BBDuk: Adapter & PhiX removal         │
│    • Trim Illumina adapters (k-mer, k=23/11)    │
│    • Remove PhiX174 spike-in reads              │
│    • Trim poly-G tails (≥5 nt)                  │
│    • Trim paired-end read overlap artifacts      │
└─────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────┐
│  Step 2 — BBDuk: Quality & length filtering     │
│    • Quality-trim 3′ ends (Phred < 0)           │
│    • Discard reads with mean quality < 3        │
│    • Discard reads with > 3 ambiguous bases     │
│    • Discard reads < 51 bp or < 33% of original │
└─────────────────────────────────────────────────┘
         │
         ▼
Filtered paired-end FASTQ + QC statistics
```

Parameters match NMDC rqcfilter2 defaults (Clum et al. 2021).

---

## Repository contents

| File | Description |
|------|-------------|
| `metagenome_qc.slurm` | SLURM job script — processes one or more sample pairs |
| `submit_all.sh` | Wrapper to submit one job per sample pair in parallel |
| `plot_qc_stats.R` | R script to parse BBDuk output and generate QC plots |

---

## Prerequisites

**BBTools (BBDuk)** — v38 or later.

On MSI, load the available module:
```bash
module load bbmap
bbduk.sh --version
```

On other clusters without a module, install via conda:
```bash
mamba create -n bbtools -c bioconda -c conda-forge bbmap -y
conda activate bbtools
```

Or download the pre-built binary:
```bash
wget https://sourceforge.net/projects/bbmap/files/latest/download -O BBTools.tar.gz
tar -xzf BBTools.tar.gz
export PATH="$PWD/bbmap:$PATH"
```

**R packages** (for visualization only):
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "forcats", "patchwork"))
```

---

## Setup

Clone and configure the SLURM script for your cluster before first use. Open `metagenome_qc.slurm` and edit the top section:

```bash
# If bbduk.sh is on PATH (module or conda), leave empty:
BBTOOLS_DIR=""

# SLURM partition name for your cluster:
#SBATCH --partition=msismall

# Notification email:
#SBATCH --mail-user=your@email.edu
```

Everything else (memory, threads, QC parameters) has sensible defaults and does not need to change for most use cases.

---

## Usage

### Single sample pair

```bash
sbatch metagenome_qc.slurm \
  --r1 /path/to/sample_R1_001.fastq.gz \
  --r2 /path/to/sample_R2_001.fastq.gz \
  --output-dir /path/to/output
```

### Whole folder (auto-paired)

The script auto-detects R1/R2 pairs by `_R1_`/`_R2_` or `_1.`/`_2.` suffixes and processes all pairs sequentially within a single job:

```bash
sbatch metagenome_qc.slurm \
  --input-dir /path/to/fastq_folder \
  --output-dir /path/to/output
```

### Multiple samples in parallel (recommended for large batches)

Use the wrapper script to submit one independent SLURM job per sample. All jobs run simultaneously subject to cluster availability:

```bash
bash submit_all.sh \
  --input-dir /path/to/fastq_folder \
  --output-dir /path/to/output
```

Monitor progress:
```bash
squeue -u $USER
```

Expects files named `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` (standard Illumina output). Edit the glob pattern in `submit_all.sh` if your naming differs.

---

## QC parameters

All parameters are set at the top of `metagenome_qc.slurm` and match NMDC rqcfilter2 defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_LEN` | 51 | Minimum read length to retain (bp) |
| `MIN_LEN_FRAC` | 0.33 | Minimum fraction of original read length to retain |
| `MIN_AVG_QUAL` | 3 | Minimum mean Phred quality score; reads below this are discarded |
| `MAX_NS` | 3 | Maximum number of ambiguous (N) bases allowed per read |
| `TRIM_POLYG` | 5 | Minimum poly-G tail length to trim |
| `JAVA_MEM` | 12g | Java heap size; BBDuk uses ~200 MB for adapter+PhiX k-mers |
| `THREADS` | 16 | CPU threads (inherited from `--cpus-per-task`) |

---

## Output structure

```
output_dir/
├── filtered/           ← cleaned read pairs for downstream analysis
│   ├── SAMPLE_R1_filtered.fastq.gz
│   └── SAMPLE_R2_filtered.fastq.gz
├── stats/              ← per-sample QC metrics
│   ├── SAMPLE_filterStats.txt    (step 1: adapter/PhiX counts per reference)
│   ├── SAMPLE_filterStats2.txt   (step 2: quality/length filter counts)
│   ├── SAMPLE_step1.log          (full BBDuk console output, step 1)
│   └── SAMPLE_step2.log          (full BBDuk console output, step 2)
├── histograms/         ← read-level distributions (post-QC)
│   ├── SAMPLE_aqhist.txt         (per-read mean quality score distribution)
│   └── SAMPLE_lhist.txt          (read length distribution)
└── logs/               ← SLURM stdout/stderr, one file per job
    ├── SAMPLE_12345.out
    └── SAMPLE_12345.err
```

The `filtered/` directory contains the files used for all downstream analyses. Everything else is QC reporting.

---

## Interpreting QC statistics

The key summary metrics appear at the bottom of each `*_step1.log` and `*_step2.log`:

**Step 1 (adapter/PhiX removal)**

| Field | What it means |
|-------|---------------|
| `KTrimmed` | Reads that had adapter sequence trimmed from the 3′ end. 20–30% is typical for short-insert Illumina libraries. Reads are retained after trimming; only counted if trimming removes enough to trigger length filters. |
| `Trimmed by overlap` | Pairs trimmed by detecting read-through into the partner read. Common for very short inserts. |
| `Total Removed` | Reads fully discarded (typically <0.1% — only reads that are adapter-only or fail length after trimming). |

**Step 2 (quality/length filtering)**

| Field | What it means |
|-------|---------------|
| `QTrimmed` | Reads that had low-quality bases trimmed from 3′ end. |
| `Low quality discards` | Reads discarded for failing mean quality or N-content thresholds. |
| `Total Removed` | Reads discarded in this step; typically 0.1–0.5% for good-quality data. |
| `Result` | Final retained reads — this is the number that proceeds to downstream analysis. |

A well-QC'd metagenome library typically retains >99% of reads, with adapter trimming occurring on 15–30% of reads (bases removed, not reads discarded).

---

## Visualizing QC results

After all jobs complete, generate summary plots across all samples:

```bash
Rscript plot_qc_stats.R /path/to/output_dir
```

Writes to `output_dir/plots/`:

| File | Description |
|------|-------------|
| `read_retention.pdf` | Stacked bar chart: retained vs. removed reads per sample |
| `removal_by_step.pdf` | Side-by-side bars: % removed in step 1 vs. step 2 |
| `length_distributions.pdf` | Post-QC read length curves, all samples overlaid |
| `quality_distributions.pdf` | Post-QC per-read mean quality curves, all samples overlaid |
| `qc_summary.pdf` | All four panels on a single page |
| `read_retention_summary.csv` | Tabular summary of read counts and retention rates |

---

## Performance

Benchmarked on MSI (Minnesota Supercomputing Institute), `msismall` partition, 16 CPUs, 16 GB RAM:

| Metric | Value |
|--------|-------|
| Input file size | ~6 GB gzipped per R1/R2 pair |
| Input read count | ~90–130 million read pairs |
| Step 1 runtime | ~4 min |
| Step 2 runtime | ~2–3 min |
| **Total per sample** | **~6–7 min** |
| Memory actually used | ~200 MB (JVM allocates 12 GB but uses ~200 MB) |

Throughput is I/O-bound (gzip compression/decompression), not CPU-bound. Increasing threads beyond 16 yields diminishing returns. For large batches, run samples in parallel rather than allocating more cores per sample.

---

## Adapting for other clusters

The script works on any SLURM cluster with BBTools available. The changes needed are:

1. **`#SBATCH --partition`** — set to your cluster's queue name
2. **`#SBATCH --mail-user`** — your email address
3. **`module load bbmap`** — replace with the appropriate module name, or remove and set `BBTOOLS_DIR` to the BBTools installation path
4. **`JAVA_MEM`** — keep below the `--mem` allocation

---

## References

**BBTools / BBDuk**
> Bushnell, B. BBTools software package. Lawrence Berkeley National Laboratory. Available at: https://sourceforge.net/projects/bbmap/

**NMDC Metagenome Reads QC Workflow** (rqcfilter2 parameters)
> Clum, A., Huntemann, M., Bushnell, B., Foster, B., Roux, S., Hajek, P. P., ... & Eloe-Fadrosh, E. A. (2021). DOE JGI Metagenome Workflow. *mSystems*, 6(2), e00804-20. https://doi.org/10.1128/mSystems.00804-20

**Shotgun metagenomics — general methods reference**
> Quince, C., Walker, A. W., Simpson, J. T., Loman, N. J., & Segata, N. (2017). Shotgun metagenomics, from sampling to analysis. *Nature Biotechnology*, 35(9), 833–844. https://doi.org/10.1038/nbt.3935

**Illumina poly-G artifact (two-color chemistry)**
> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
