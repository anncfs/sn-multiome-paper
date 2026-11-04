# Additional UMAP plots for SI

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PARAMS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
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

# config
source(args$config)

# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

# output dir
dir_out <- file.path(dir_scripts, "data/figures-supplementary/umaps")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# sc data
sc <- readRDS(file.path(dir_manuscript, "sc_rna.rds"))
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))


# -----------------------------------------------------------------------
# FigSX: brain region UMAP
# -----------------------------------------------------------------------

# UMAP coloured by brain region
fig_umap_region <- DimPlot(sc, cols = named_colours$region, group.by = "region",
                      reduction = "umap.harmony", raster = FALSE,
                      label = FALSE, label.size = 6, pt.size = 0.1) +
  scale_color_manual(values = named_colours$region,
                     labels = gsub("PFC", "Prefrontal cortex", names(named_colours$region))) +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())

# fig_umap_region
ggsave(file.path(dir_out, "fig_umap_region.pdf"), fig_umap_region, width = 8, height = 6)


# -----------------------------------------------------------------------
# FigSX: UMAP WNN
# -----------------------------------------------------------------------

# umaps based on ATAC or WNN graphs
umap_list <- sapply(c("umap.harmony.atac", "umap.wnn"), function(i) {
  fig_umap <- DimPlot(sc, cols = named_colours$celltype_name, group.by = "celltype_name",
                         reduction = i, raster = FALSE,
                         label = FALSE, label.size = 6, pt.size = 0.1) +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())
  return(fig_umap)
})

# relative RNA weight per cell for WNN
fig_weighting <- data.frame(RNA.weight = sc$RNA.weight) |>
  arrange(desc(RNA.weight)) |>
  mutate(idx = row_number()) |>
  ggplot(aes(x = idx, y = RNA.weight)) +
    geom_point(size = 0.3, color = "black") +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.margin = margin(t = 30, r = 30, b = 30, l = 30)) +
  xlab("Cell index") +
  ylab("RNA weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(0, 150000))

# celltype legend
fig_legend <- data.frame(x = 1, y = seq_along(named_colours$celltype_name),
                         celltype = factor(names(named_colours$celltype_name),
                                           levels = names(named_colours$celltype_name))) |>
  ggplot(aes(x = x, y = y, color = celltype)) +
    geom_point(size = 4) + theme_void() +
    scale_color_manual(values = named_colours$celltype_name,
                     labels = names(named_colours$celltype_name)) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_blank())

ct_legend <- wrap_elements(full = cowplot::get_legend(fig_legend))

design <- "
  123
  124"
umap_list[[1]] + umap_list[[2]] + ct_legend + fig_weighting +
  plot_layout(design = design, widths = c(1,1,0.4), heights = c(0.6, 0.4))

ggsave(file.path(dir_out, "fig_map_wnn.pdf"), width = 13, height = 6)

# check percentiles where RNA weight more important
df_wnn <- sc@meta.data |>
  select(RNA.weight, AMPFRAG.weight)

cat(
  sprintf("%.1f %% cells have RNA.weight >= 0.5", 100 * sum(df_wnn$RNA.weight >= 0.5) / nrow(df_wnn)),
  file = file.path(dir_out, "wnn-stats.txt"))
