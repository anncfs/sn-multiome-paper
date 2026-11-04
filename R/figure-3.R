# Figure 3 plots (heritability)

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PARAMS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(argparse)
})

# parse arguments
parser <- ArgumentParser(description = "snMultiome figures")
parser$add_argument(
  "--config", required = TRUE, type = "character")

args <- parser$parse_args()

# config for file paths
source(args$config)

# functions
source(file.path(dir_scripts, "analysis/analysis-seurat/R/functions-figures.R"))


# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

# output dir
dir_out <- file.path(dir_scripts, "data/figures-main/figure-3")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# combined heritability enrichment results
df_heritability <- readRDS(file.path(dir_manuscript, "heritability_enrichment_results.rds")) |>
  mutate(CELLTYPE = gsub("-", "_", CELLTYPE))

# colours
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# GWAS info
df_gwas <- read.csv(file.path(dir_manuscript, "gwas_info.csv"), header = T) |>
  as_tibble() |>
  filter(Plot == TRUE) |>
  mutate(Category = factor(Category, levels = names(named_colours$Category))) |>
  arrange(Category, GWAS)
df_gwas$GWAS <- factor(df_gwas$GWAS, levels = df_gwas$GWAS)

# -----------------------------------------------------------------------
# A: RNA MAGMA SUBCELLTYPE
# -----------------------------------------------------------------------

# unconditioned subcelltype analysis
df <- df_heritability |>
  filter(ANALYSIS == "unconditioned", GROUPING == "subcelltype") |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$subcelltype))) |>
  filter(GWAS %in% filter(df_gwas, Plot == TRUE)$GWAS)

fig3A <- plot_heatmap_enrichment(df, df_gwas = df_gwas)
# fig3A

pdf(file.path(dir_out, "fig3A_magma.pdf"), width = 13, height = 10)
fig3A
invisible(dev.off())



# -----------------------------------------------------------------------
# C: ATAC LDSC CELLTYPE
# -----------------------------------------------------------------------

# LDSC ATAC analysis
df <- df_heritability |>
  filter(METHOD == "LDSC", GROUPING == "celltype", MODALITY == "ATAC") |>
  mutate(CELLTYPE = factor(CELLTYPE)) |>
  filter(GWAS %in% filter(df_gwas, Plot == TRUE)$GWAS)

# udpate celltype names and order
levels(df$CELLTYPE) <- sapply(seq_along(levels(df$CELLTYPE)), function(x) {
  grep(levels(df$CELLTYPE)[x], names(named_colours$celltype_name), value = T)
})

df <- df |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$celltype_name)))
  
fig3C <- plot_heatmap_enrichment(df,
                                 coefficient_name = "TAU.NORM", 
                                 df_gwas = df_gwas,
                                 colour = "royalblue4")

pdf(file.path(dir_out, "fig3C_atac.pdf"), width = 13, height = 6)
fig3C
dev.off()

pdf(file.path(dir_out, "fig3AC_combined.pdf"), width = 13, height = 14)
fig3A@bottom_annotation <- NULL
draw(fig3A %v% fig3C, merge_legend = FALSE, ht_gap = unit(1.5, "cm"), 
     legend_grouping = "original")
invisible(dev.off())


# -----------------------------------------------------------------------
# B: CONDITIONAL ANALYSES
# -----------------------------------------------------------------------

gwas_to_plot <- c("BMI_2018", "BP_2025", "MD_2025", "SCZ_2022", "AAM_2024")
names(gwas_to_plot) <- c("BMI", "Bipolar\ndisorder", "Major\ndepression", "Schizophrenia", "Age at\nmenarche")

# select iterative analyses for GWAS of interest
df <- df_heritability |>
  filter(ANALYSIS == "iterative", GWAS %in% gwas_to_plot) |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$subcelltype))) |>
  select(CELLTYPE, GWAS, FDR, BETA, ROUND, NEURON)

max_beta <- ceiling(max(abs(df$BETA)))

plist <- lapply(seq_along(gwas_to_plot), function(i) {
  p <- plot_heatmap_conditional_iterative(df, id = gwas_to_plot[i], max_beta = max_beta, 
                                          star = TRUE, title = names(gwas_to_plot)[i])
  if (i > 1) {
    p@heatmap_param$show_heatmap_legend <- FALSE
  }
  return(p)
})

# add custom legends
# circle for conditioning variables
lgd_list <- list(
  Legend(labels = "", title = "Conditioning\nvariable", type = "points", pch = 16,
         title_gp = gpar(fontsize = 12, fontface = "bold"),
         size = unit(0.4, "cm"), legend_gp = gpar(col = "blue"), background = "white"),
  Legend(title = "FDR", labels = c("< 0.05", "< 0.01"),
         title_gp = gpar(fontsize = 12, fontface = "bold"),
         labels_gp = gpar(fontsize = 12),
         graphics = list(function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 14)),
                         function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 14)))))

pdf(paste0(dir_out, "/fig3B_iterative.pdf"), width = 12, height = 10)
to_plot <- Reduce("+", plist)
fig3B <- draw(to_plot, ht_gap = unit(0.5, "cm"), merge_legend = TRUE, annotation_legend_list = lgd_list)
invisible(dev.off())

