# Supplementary Figures/Tables for Figure 1 (ATAC peaks)

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PARAMS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
suppressPackageStartupMessages({
  # library(Seurat)
  # library(Signac)
  # library(EnsDb.Hsapiens.v86)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(UpSetR)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(argparse)
})

# parse arguments
parser <- ArgumentParser(description = "snMultiome figures")
parser$add_argument(
  "--config", required = TRUE, 
  help = "Path to config.R", 
  type = "character")

args <- parser$parse_args()

source(args$config)

# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

dir_in <- file.path(dir_data, "15_manuscript_figures/data")
dir_out <- file.path(dir_scripts, "data/figures-supplementary/atac-peaks")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# load peak set
peaks <- readRDS(file.path(dir_in, "peakset-AMPFRAG.rds"))

# colours
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))


# -----------------------------------------------------------------------
# FIGURE SX: PEAK CELL TYPE SPECIFICITY
# -----------------------------------------------------------------------

df_peaks <- as.data.frame(peaks) |>
  mutate(peak_id = paste(seqnames, start, end, sep = "-")) |>
  select(peak_id, peak_called_in)

# matrix for upset plot / cell types per peak
df_long <- df_peaks |>
  tidyr::separate_rows(peak_called_in, sep = ",") |>
  mutate(celltype = trimws(peak_called_in))

ct_matrix <- df_long |>
  mutate(value = 1L) |>
  tidyr::pivot_wider(
    id_cols = peak_id, 
    names_from = celltype,
    values_from = value,
    values_fill = 0L)

# set cell type names and order
ct_names <- data.frame(orig = colnames(ct_matrix)[-1]) |>
  mutate(new = tools::toTitleCase(gsub("-", " ", orig))) |>
  mutate(new = recode(new, "Oligodendrocyte Precursor" = "OPC"))
colnames(ct_matrix) <- c("peak_id", ct_names$new)

ct_order <- names(named_colours$celltype_name)

# Upset plot of peak overlaps
p_input <- ct_matrix |>
  select(all_of(ct_order)) |>
  as.data.frame()

p <- upset(p_input,
           sets = rev(ct_order),
           keep.order = TRUE,
           nsets = length(ct_order),
           nintersects = 20,
           order.by = "freq",
           sets.bar.color = rev(named_colours$celltype_name),
           main.bar.color = "gray40",
           mainbar.y.label = "Number of peaks",
           sets.x.label = "Peaks per cell type",
           matrix.color = "gray40",
           number.angles = 0,
           text.scale = c(1.8,1.5,1.5,1.5,1.5,1.1),
           mb.ratio = c(0.6,0.4))

pdf(file.path(dir_out, "fig_upset_peak_counts.pdf"), width = 10, height = 6)
print(p)
invisible(dev.off())


# -----------------------------------------------------------------------
# COUNTS TABLE: CELL TYPE SPECIFICITY
# -----------------------------------------------------------------------

peak_counts <- df_peaks |>
  summarise(n = n(), .by = peak_called_in) |>
  arrange(desc(n))
n_total <- sum(peak_counts$n)

# peaks specific to 1 cell type
counts_1x <- peak_counts |>
  filter(peak_called_in %in% ct_names$orig)
n_1x <- sum(counts_1x$n)

# neuronal specific peaks
counts_neu <- peak_counts |>
  filter(peak_called_in %in% c("excitatory-neuron", "inhibitory-neuron",
                               "excitatory-neuron,inhibitory-neuron",
                               "inhibitory-neuron,excitatory-neuron"))
n_neu <- sum(counts_neu$n)

# save peak counts 
writeLines(c(
  sprintf("%-20s %d ", "Total peaks:", n_total),
  sprintf("%-20s %d (%.1f %%)", "1 cell type peaks:", n_1x, 100*(n_1x/n_total)),
  sprintf("%-20s %d (%.1f %%)", "Neuronal peaks:", n_neu, 100*(n_neu/n_total))
), con = file.path(dir_out, "peak_counts_summary.txt"))
