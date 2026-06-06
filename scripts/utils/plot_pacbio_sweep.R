#!/usr/bin/env Rscript
# plot_pacbio_sweep.R
# Generates ROC curve and two stacked bar charts from PacBio sweep outputs.
#
# Usage:
#   Rscript plot_pacbio_sweep.R <sweep_results.tsv> <family_counts.tsv> <output_dir>
#
# Outputs:
#   <output_dir>/roc.png
#   <output_dir>/bar_rrna.png
#   <output_dir>/bar_nonrrna.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_pacbio_sweep.R <sweep_results.tsv> <family_counts.tsv> <output_dir>")
}

sweep_tsv  <- args[1]
family_tsv <- args[2]
out_dir    <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- ROC curve ---
sweep <- read.delim(sweep_tsv, stringsAsFactors = FALSE) %>%
  mutate(selectivity = 1 - fpr) %>%
  arrange(num_seeds)

p_roc <- ggplot(sweep, aes(x = selectivity, y = sensitivity, label = num_seeds)) +
  geom_line(colour = "grey60", linetype = "dashed") +
  geom_point(aes(colour = num_seeds), size = 2.5) +
  geom_text(nudge_x = 0.004, nudge_y = 0.004, size = 2.5, colour = "grey30") +
  scale_colour_viridis_c(name = "num_seeds", direction = -1) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "forestgreen") +
  scale_x_continuous(limits = c(0, 1.02)) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(x = "Selectivity (1 - FPR)", y = "Sensitivity",
       title = "SortMeRNA PacBio sweep: ROC",
       subtitle = "Labels = num_seeds value") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "roc.png"), p_roc, width = 7, height = 5.5, dpi = 300)
message("  Saved roc.png")

# --- Stacked bar charts ---
fam <- read.delim(family_tsv, stringsAsFactors = FALSE) %>%
  mutate(subunit = case_when(
    grepl("SSU.*Bacteria|SSU.*Archaea", family)   ~ "16S",
    grepl("SSU.*Eukaryota",             family)   ~ "18S",
    grepl("LSU.*Bacteria|LSU.*Archaea", family)   ~ "23S",
    grepl("LSU.*Eukaryota",             family)   ~ "28S",
    grepl("5\\.8S",                     family)   ~ "5.8S",
    grepl("5S",                         family)   ~ "5S",
    TRUE                                          ~ "Unknown"
  )) %>%
  group_by(num_seeds, rna_type, subunit) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    subunit   = factor(subunit, levels = c("16S", "18S", "23S", "28S", "5S", "5.8S", "Unknown", "No alignment")),
    num_seeds = factor(num_seeds, levels = sort(unique(num_seeds)))
  )

palette <- c(
  "16S"     = "#1f77b4",
  "18S"     = "#17becf",
  "23S"     = "#d62728",
  "28S"     = "#ff7f0e",
  "5S"      = "#2ca02c",
  "5.8S"    = "#98df8a",
  "Unknown"      = "#7f7f7f",
  "No alignment" = "#ffffff"
)

make_bar <- function(data, type_label, title) {
  ggplot(data %>% filter(rna_type == type_label),
         aes(x = num_seeds, y = count, fill = subunit)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette, name = "rRNA subunit", drop = FALSE) +
    labs(x = "num_seeds", y = "Aligned reads", title = title) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_rrna <- make_bar(fam, "rrna",    "rRNA reads: aligned by subunit vs num_seeds")
p_non  <- make_bar(fam, "nonrrna", "Non-rRNA reads (FP): aligned by subunit vs num_seeds")

ggsave(file.path(out_dir, "bar_rrna.png"),    p_rrna, width = 9, height = 5, dpi = 300)
ggsave(file.path(out_dir, "bar_nonrrna.png"), p_non,  width = 9, height = 5, dpi = 300)
message("  Saved bar_rrna.png and bar_nonrrna.png")
