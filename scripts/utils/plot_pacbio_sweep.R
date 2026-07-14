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
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_pacbio_sweep.R <sweep_results.tsv> <family_counts.tsv> <output_dir>")
}

sweep_tsv  <- args[1]
family_tsv <- args[2]
out_dir    <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- ROC curve (3 panels, one per e-value, auto-scaled axes) ---
sweep <- read.delim(sweep_tsv, stringsAsFactors = FALSE) %>%
  mutate(
    selectivity = 1 - fpr,
    label       = paste0("ms=", min_lis, ",ns=", num_seeds),
    evalue      = factor(evalue, levels = sort(unique(evalue)))
  )

# Use the most permissive e-value panel (largest numeric e-value = widest x-range)
# for the shared x-axis limits so all panels are comparable.
max_ev <- max(as.numeric(levels(sweep$evalue)))
x_limits <- sweep %>%
  filter(abs(as.numeric(as.character(evalue)) - max_ev) < max_ev * 1e-6) %>%
  summarise(lo = min(selectivity) - 0.02, hi = 1.0)

p_roc <- ggplot(sweep, aes(x = selectivity, y = sensitivity)) +
  geom_point(size = 2, colour = "#2c3e50") +
  geom_label_repel(
    data = sweep %>% filter(abs(as.numeric(as.character(evalue)) - max_ev) < max_ev * 1e-6),
    aes(label = label), size = 2.2, colour = "grey20",
    box.padding = 0.3, point.padding = 0.2,
    max.overlaps = Inf, min.segment.length = 0) +
  facet_wrap(~ evalue, scales = "free_y", ncol = 3,
             labeller = labeller(evalue = function(x) paste0("e = ", x))) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "forestgreen") +
  scale_x_continuous(limits = c(x_limits$lo, x_limits$hi)) +
  labs(x = "Selectivity (1 - FPR)", y = "Sensitivity",
       title = "SortMeRNA PacBio sweep: ROC by e-value",
       subtitle = "Labels: ms = min_lis, ns = num_seeds") +
  theme_bw(base_size = 11) +
  theme(panel.spacing = unit(1, "lines"))

ggsave(file.path(out_dir, "roc.png"), p_roc, width = 14, height = 5, dpi = 300)
message("  Saved roc.png")

# --- Stacked bar charts ---
fam_raw <- read.delim(family_tsv, stringsAsFactors = FALSE) %>%
  mutate(subunit = case_when(
    grepl("SSU.*Bacteria|SSU.*Archaea", family, ignore.case = TRUE)   ~ "16S",
    grepl("SSU.*Eukaryota",             family, ignore.case = TRUE)   ~ "18S",
    grepl("LSU.*Bacteria|LSU.*Archaea", family, ignore.case = TRUE)   ~ "23S",
    grepl("LSU.*Eukaryota",             family, ignore.case = TRUE)   ~ "28S",
    grepl("5.8S",                       family, ignore.case = TRUE)   ~ "5.8S",
    grepl("5S",                         family, ignore.case = TRUE)   ~ "5S",
    TRUE                                                              ~ "Unknown"
  ))

palette <- c(
  "16S"          = "#1f77b4",
  "18S"          = "#17becf",
  "23S"          = "#d62728",
  "28S"          = "#ff7f0e",
  "5S"           = "#2ca02c",
  "5.8S"         = "#98df8a",
  "Unknown"      = "#7f7f7f",
  "No alignment" = "#ffffff"
)

subunit_levels <- c("16S", "18S", "23S", "28S", "5S", "5.8S", "Unknown", "No alignment")

make_bar <- function(data, x_var, type_label, title, xlab) {
  df <- data %>%
    filter(rna_type == type_label) %>%
    group_by(across(all_of(c(x_var, "subunit")))) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(
      subunit      = factor(subunit, levels = subunit_levels),
      x_val        = factor(.data[[x_var]], levels = sort(unique(.data[[x_var]])))
    )
  ggplot(df, aes(x = x_val, y = count, fill = subunit)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette, name = "rRNA subunit", drop = FALSE) +
    labs(x = xlab, y = "Aligned reads", title = title) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

make_row <- function(x_var, xlab, row_title, fix_ev = NULL, fix_ns = NULL, fix_lis = NULL) {
  filtered <- fam_raw
  if (!is.null(fix_ev))  filtered <- filtered %>% filter(abs(evalue_num - fix_ev) < fix_ev * 1e-6)
  if (!is.null(fix_ns))  filtered <- filtered %>% filter(num_seeds == fix_ns)
  if (!is.null(fix_lis)) filtered <- filtered %>% filter(min_lis   == fix_lis)
  p1 <- make_bar(filtered, x_var, "rrna",    paste0(row_title, " - rRNA"),    xlab)
  p2 <- make_bar(filtered, x_var, "nonrrna", paste0(row_title, " - non-rRNA"), xlab)
  p1 + p2
}

# Fix reference values for non-swept parameters
fam_raw <- fam_raw %>% mutate(evalue_num = as.numeric(evalue))
ref_ev_num <- max(fam_raw$evalue_num)
ref_ns  <- "2"   # num_seeds reference
ref_lis <- "2"   # min_lis reference

row_ns  <- make_row("num_seeds", "num_seeds",
                    paste0("vs num_seeds  [evalue=", ref_ev_num, ", min_lis=", ref_lis, "]"),
                    fix_ev = ref_ev_num, fix_lis = as.integer(ref_lis))
row_lis <- make_row("min_lis",   "min_lis",
                    paste0("vs min_lis  [evalue=", ref_ev_num, ", num_seeds=", ref_ns, "]"),
                    fix_ev = ref_ev_num, fix_ns = as.integer(ref_ns))
row_ev  <- make_row("evalue_num", "e-value",
                    paste0("vs e-value  [num_seeds=", ref_ns, ", min_lis=", ref_lis, "]"),
                    fix_ns = as.integer(ref_ns), fix_lis = as.integer(ref_lis))

combined <- (row_ns / row_lis / row_ev) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "rRNA family breakdown (left: rRNA reads, right: non-rRNA reads)") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "bar_combined.png"), combined, width = 14, height = 14, dpi = 300)
message("  Saved bar_combined.png")
