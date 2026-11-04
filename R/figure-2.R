# Figure 2 plots

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
dir_out <- file.path(dir_scripts, "data/figures-main/figure-2")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# sc data
sc <- readRDS(file.path(dir_manuscript, "sc_rna.rds"))
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# variance partition results
vp_stack <- readRDS(file.path(dir_manuscript, "vp_dreamlet_stacked.rds"))
vp_celltype <- readRDS(file.path(dir_manuscript, "vp_dreamlet_celltype.rds"))


# -----------------------------------------------------------------------
# VARIANCE PARTITION RESULTS
# -----------------------------------------------------------------------

vp_stack <- vp_stack |>
  select(!assay) |>
  tidyr::pivot_longer(cols = !gene, names_to = "covariate", values_to = "variance") |>
  mutate(covariate = replace_values(stringr::str_to_title(covariate),
                                      "Celltype" ~ "Cell type",
                                      "Scale.pmd" ~ "PMD (scaled)",
                                      "Scale.age" ~ "Age (scaled)",
                                      "Brainbank" ~ "Brain bank",
                                      "Log.nfeature_rna" ~ "No. of genes (log)",
                                      "Percent.mt" ~ "% mitochondrial",
                                      "Percent.ribosomal" ~ "% ribosomal"))

# variance explained by covariate
vp_mean_order <- vp_stack |>
  group_by(covariate) |>
  summarise(mean = mean(variance)) |>
  arrange(desc(mean)) |>
  dplyr::filter(covariate != "Residuals") |>
  pull(covariate)

# vp (stacked) results
fig2_vp <- vp_stack |>
  mutate(covariate = factor(covariate, levels = c(vp_mean_order, "Residuals"))) |>
  ggplot(aes(x = covariate, y = 100 * variance)) +
    geom_violin(scale = "width", fill = "grey80", color = NA) +
    theme_classic() +
    ylab("Total variance explained (%)") +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.colour = "black", outlier.size = 0.5) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          legend.position = "none",
          aspect.ratio = 0.8,
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12), 
          axis.title.y = element_text(size = 14)) 

# fig2_vp
ggsave(file.path(dir_out, "fig2_varpart.pdf"), fig2_vp, width = 5, height = 5)

# vp (for region) by cell type
vp_celltype <- vp_celltype |>
  mutate(celltype = factor(assay)) 
levels(vp_celltype$celltype) <- data.frame(celltype = levels(vp_celltype$celltype)) |>
  mutate(celltype = tools::toTitleCase(gsub("-", " ", celltype))) |>
  mutate(celltype = recode(celltype, "Oligodendrocyte Precursor" = "OPC")) |>
  pull(celltype)
vp_celltype <- vp_celltype |>
  mutate(celltype = factor(celltype, levels = names(named_colours$celltype_name)))

fig2_region <- vp_celltype |>
  ggplot(aes(x = celltype, y = 100 * region)) +
  geom_violin(scale = "width", fill = "grey80", color = NA) +
  theme_classic() +
  ylab("Variance explained by region (%)") +
  geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.colour = "black", outlier.size = 0.5) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14)) 

# fig2_region

ggsave(file.path(dir_out, "fig2_region.pdf"), fig2_region, width = 5, height = 5)


# -----------------------------------------------------------------------
# UMAP - SUBCELLTYPE
# -----------------------------------------------------------------------

# UMAP coloured by sub cell type
fig2_umap <- DimPlot(sc, cols = named_colours$subcelltype, group.by = "subcelltype",
                      reduction = "umap.harmony", raster = FALSE,
                      label = FALSE, pt.size = 0.1) +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())

# fig2_umap
ggsave(file.path(dir_out, "fig2_umap.pdf"), fig2_umap, width = 6, height = 6)


# -----------------------------------------------------------------------
# PROPORTIONS - SUBCELLTYPE
# -----------------------------------------------------------------------

# select metadata columns
df_annot <- sc@meta.data |>
  mutate(Region = recode(region, "PFC" = "Prefrontal cortex")) |>
  select(subcelltype, cellclass, Region)
named_colours$Region <- setNames(named_colours$region,
                                 nm = gsub("PFC", "Prefrontal cortex", names(named_colours$region)))

fig2_proportions <- plot_barplot_vertical(df_annot, named_colours,
                                          "subcelltype", "Region",
                                          group_colours = T) + 
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)")
  


ggsave(file.path(dir_out, "fig2_proportions.pdf"), fig2_proportions, width = 6, height = 6)



