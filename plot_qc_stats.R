library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)
library(patchwork)

# =============================================================================
# Metagenome QC — Summary visualization
#
# Reads BBDuk output from the stats/ and histograms/ subdirectories produced
# by metagenome_qc.slurm and generates:
#   1. Per-sample read retention summary (step 1 + step 2)
#   2. Read length distributions (post-QC)
#   3. Average quality score distributions (post-QC)
#
# Usage:
#   Rscript plot_qc_stats.R /path/to/output_dir
#   # or source interactively and set OUT_DIR below
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else "."

STATS_DIR <- file.path(OUT_DIR, "stats")
HIST_DIR  <- file.path(OUT_DIR, "histograms")
PLOT_DIR  <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

# =============================================================================
# 1.  Parse step logs → per-sample read counts
# =============================================================================

parse_bbduk_log <- function(path) {
  lines <- readLines(path, warn = FALSE)

  grab <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE)
    if (length(hit) == 0) return(NA_real_)
    as.numeric(sub(".*?(\\d+) reads.*", "\\1", hit[1]))
  }

  list(
    input         = grab("^Input:"),
    ktrimmed      = grab("^KTrimmed:"),
    poly_trimmed  = grab("^Polymer-trimmed:"),
    overlap_trim  = grab("^Trimmed by overlap:"),
    total_removed = grab("^Total Removed:"),
    result        = grab("^Result:")
  )
}

step1_logs <- list.files(STATS_DIR, pattern = "_step1\\.log$", full.names = TRUE)
step2_logs <- list.files(STATS_DIR, pattern = "_step2\\.log$", full.names = TRUE)

if (length(step1_logs) == 0) stop("No step1 log files found in: ", STATS_DIR)

parse_logs <- function(logs, step_label) {
  map_dfr(logs, function(f) {
    sample <- basename(f) |> sub(pattern = "_step[12]\\.log$", replacement = "")
    m <- parse_bbduk_log(f)
    tibble(
      sample        = sample,
      step          = step_label,
      input         = m$input,
      ktrimmed      = m$ktrimmed,
      poly_trimmed  = m$poly_trimmed,
      overlap_trim  = m$overlap_trim,
      total_removed = m$total_removed,
      result        = m$result,
      pct_retained  = result / input * 100
    )
  })
}

s1 <- parse_logs(step1_logs, "Step 1: Adapter/PhiX")
s2 <- parse_logs(step2_logs, "Step 2: Quality filter")

# Combined retention: reads surviving both steps vs original input
combined <- s1 |>
  select(sample, input_reads = input) |>
  left_join(s2 |> select(sample, final_reads = result), by = "sample") |>
  mutate(
    reads_removed = input_reads - final_reads,
    pct_retained  = final_reads / input_reads * 100
  ) |>
  arrange(sample)

cat("\n=== Read retention summary ===\n")
combined |>
  mutate(across(where(is.numeric), \(x) formatC(x, format = "f", digits = 1, big.mark = ","))) |>
  print(n = Inf)

# =============================================================================
# 2.  Plot: stacked bar — reads retained vs removed, per sample
# =============================================================================

plot_data <- combined |>
  pivot_longer(c(final_reads, reads_removed),
               names_to = "category", values_to = "reads") |>
  mutate(
    category = factor(category,
                      levels = c("reads_removed", "final_reads"),
                      labels = c("Removed", "Retained")),
    reads_M = reads / 1e6,
    sample  = fct_reorder(sample, -reads_M)
  )

label_y_offset <- max(combined$input_reads) / 1e6 * 0.02
combined_labeled <- combined |> mutate(label_y = input_reads / 1e6 + label_y_offset)

p_bar <- ggplot(plot_data, aes(x = sample, y = reads_M, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(
    data = combined_labeled,
    aes(x = sample, y = label_y, label = sprintf("%.1f%%", pct_retained), fill = NULL),
    size = 2.8, hjust = 0.5
  ) +
  scale_fill_manual(values = c("Retained" = "#2166ac", "Removed" = "#d73027")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Read retention after QC",
    subtitle = "Adapter/PhiX removal + quality filtering (BBDuk)",
    x = NULL, y = "Reads (millions)", fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank())

# Step-by-step removal breakdown
step_data <- bind_rows(s1, s2) |>
  mutate(
    pct_removed = total_removed / input * 100,
    step = factor(step, levels = c("Step 1: Adapter/PhiX", "Step 2: Quality filter")),
    sample = fct_reorder(sample, -pct_removed)
  )

p_steps <- ggplot(step_data, aes(x = sample, y = pct_removed, fill = step)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Step 1: Adapter/PhiX"  = "#4dac26",
                               "Step 2: Quality filter" = "#d01c8b")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Reads removed per QC step (%)",
    x = NULL, y = "% of input reads removed", fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank())

# =============================================================================
# 3.  Parse histograms
# =============================================================================

parse_hist <- function(dir, pattern, col_names) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)

  map_dfr(files, function(f) {
    sample <- basename(f) |>
      sub(pattern = pattern, replacement = "") |>
      sub(pattern = "\\.txt$", replacement = "")
    d <- read_tsv(f, comment = "#", col_names = col_names,
                  col_types = cols(.default = "d"), show_col_types = FALSE)
    # lhist has no fraction column — compute it
    if (!"fraction" %in% names(d))
      d <- mutate(d, fraction = count / sum(count))
    d$sample <- sample
    d
  })
}

# lhist: 2 columns (length, count); aqhist: 3 columns (quality, count, fraction)
lhist  <- parse_hist(HIST_DIR, "_lhist",  c("length",  "count"))
aqhist <- parse_hist(HIST_DIR, "_aqhist", c("quality", "count", "fraction"))

# Read length distribution
p_len <- NULL
if (!is.null(lhist) && nrow(lhist) > 0) {
  p_len <- lhist |>
    filter(count > 0) |>
    ggplot(aes(x = length, y = fraction * 100, colour = sample, group = sample)) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    labs(
      title = "Read length distribution (post-QC)",
      x = "Read length (bp)", y = "% of reads", colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.text = element_text(size = 7))
}

# Average quality distribution
p_qual <- NULL
if (!is.null(aqhist) && nrow(aqhist) > 0) {
  p_qual <- aqhist |>
    filter(count > 0) |>
    ggplot(aes(x = quality, y = fraction * 100, colour = sample, group = sample)) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    labs(
      title = "Per-read average quality distribution (post-QC)",
      x = "Mean Phred quality score", y = "% of reads", colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.text = element_text(size = 7))
}

# =============================================================================
# 4.  Save plots
# =============================================================================

ggsave(file.path(PLOT_DIR, "read_retention.pdf"), p_bar,
       width = max(8, nrow(combined) * 0.4 + 3), height = 5)

ggsave(file.path(PLOT_DIR, "removal_by_step.pdf"), p_steps,
       width = max(8, nrow(combined) * 0.4 + 3), height = 5)

if (!is.null(p_len))
  ggsave(file.path(PLOT_DIR, "length_distributions.pdf"), p_len, width = 10, height = 5)

if (!is.null(p_qual))
  ggsave(file.path(PLOT_DIR, "quality_distributions.pdf"), p_qual, width = 10, height = 5)

# Combined summary page
if (!is.null(p_len) && !is.null(p_qual)) {
  combined_plot <- (p_bar | p_steps) / (p_len | p_qual) +
    plot_annotation(title = "Metagenome QC Summary", theme = theme(plot.title = element_text(size = 14)))
  ggsave(file.path(PLOT_DIR, "qc_summary.pdf"), combined_plot, width = 16, height = 10)
}

# Write summary table as CSV
write_csv(combined, file.path(PLOT_DIR, "read_retention_summary.csv"))
cat("\nPlots written to:", PLOT_DIR, "\n")
