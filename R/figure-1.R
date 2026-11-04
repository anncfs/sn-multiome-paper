# Figure 1 plots

# Usage:
# Rscript figure-1.r --config /path/to/config.R

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

# config for file paths
source(args$config)

# functions
source(file.path(dir_scripts, "analysis/analysis-seurat/R/functions-figures.R"))


# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

# output dir
dir_out <- file.path(dir_scripts, "data/figures-main/figure-1")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# sc data
sc <- readRDS(file.path(dir_manuscript, "sc_rna.rds"))
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# geneactivity
mean_geneactivity <- readRDS(file.path(dir_manuscript, "mean_geneactivity_celltype.rds"))
mean_expression <- readRDS(file.path(dir_manuscript, "mean_rna_celltype.rds"))

# marker genes
marker_genes <- read.csv(file.path(dir_scripts, "data/metadata/figure-marker-genes.csv"), 
                         stringsAsFactors = F, col.names = c("celltype", "gene", "abbrev"))


# -----------------------------------------------------------------------
# Fig1A: cell type UMAP
# -----------------------------------------------------------------------

# UMAP coloured by major cell type
fig1A_umap <- DimPlot(sc, cols = named_colours$celltype_name, group.by = "celltype_name",
                 reduction = "umap.harmony", raster = FALSE,
                 label = TRUE, label.size = 6, pt.size = 0.1) +
  NoLegend() +
  theme(aspect.ratio = 1, plot.title = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())

# fig1A_umap
ggsave(file.path(dir_out, "fig1A_umap.pdf"), fig1A_umap, width = 6, height = 6)

# No. cells per celltype
cell_counts <- sc@meta.data |>
  summarise(.by = celltype_name, n = n()) |>
  arrange(desc(n)) |>
  mutate(celltype_name = factor(celltype_name, levels = unique(celltype_name))) 

fig1A_inset <- cell_counts |>
  ggplot(aes(x = n, y = celltype_name, fill = celltype_name)) +
    geom_col() +
    theme_void() +
    scale_fill_manual(values = named_colours$celltype_name) +
    geom_text(aes(label = scales::comma(n)), hjust = -0.2, size = 5) +
    NoLegend() +
    xlim(c(0, 80000)) +
    xlab("# cells") +
    theme(axis.title.x = element_text(size = 20, hjust = 0.2)) +
    theme(plot.margin = margin(0,0,20,0))
  
# fig1A_inset
ggsave(file.path(dir_out, "fig1A_inset.pdf"), fig1A_inset, width = 4, height = 4)

fig1A <- fig1A_umap + inset_element(fig1A_inset, 0, 0, 0.28, 0.28, align_to = "full")
ggsave(file.path(dir_out, "fig1A.pdf"), fig1A, width = 8, height = 8)



# -----------------------------------------------------------------------
# CELL TYPE X REGION BARPLOTS
# -----------------------------------------------------------------------

df_annot <- sc@meta.data |>
  select(celltype_name, region) |>
  mutate(celltype_name = factor(celltype_name, levels = rev(cell_counts$celltype_name)))

# # frequency
# fig1C_counts <- plot_barplot_vertical(df_annot, named_colours, "region", "celltype_name", proportion = F) +
#   scale_fill_manual(NULL, values = named_colours$celltype_name,
#                     breaks = names(named_colours$celltype_name)) +
#   scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex")) +
#   scale_x_continuous(breaks = c(0, 20000, 40000), labels = c(0, 20, 40)) +
#   xlab(bquote("Number of nuclei "(10^3)))

# proportion
fig1C_proportion <- plot_barplot_vertical(df_annot, named_colours, "region", "celltype_name", proportion = T) +
  scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex")) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name))  +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)")

# fig1C_counts + theme(legend.key.size = unit(0.45, 'cm'))
# ggsave(file.path(dir_out, "fig1C_counts.pdf"), width = 6, height = 2.5)
fig1C_proportion + theme(legend.key.size = unit(0.45, 'cm'))
ggsave(file.path(dir_out, "fig1C_proportion.pdf"), width = 6, height = 2.5)


# -----------------------------------------------------------------------
# MARKER GENE PROFILES
# -----------------------------------------------------------------------

# set plot order
gene_order <- marker_genes$gene
celltype_order <- marker_genes |> 
  distinct(celltype, abbrev)

# scaled data matrix
df1 <- LayerData(mean_expression, layer = "scale.data", features = gene_order)
df2 <- LayerData(mean_geneactivity, layer = "scale.data", features = gene_order)

df1 <- t(df1[gene_order, celltype_order$abbrev])
df2 <- t(df2[gene_order, celltype_order$abbrev])

rownames(df1) <- celltype_order$celltype
rownames(df2) <- celltype_order$celltype

h1 <- get_heatmap_expression(df1, title = "Mean gene expression (RNA)")
h2 <- get_heatmap_expression(df2, title = "Predicted gene activity (ATAC)")
h1@heatmap_param$show_heatmap_legend <- FALSE
h_list <- h1 %v% h2

fig1B <- draw(h_list)
pdf(paste0(dir_out, "/fig1B.pdf"), width = 10, height = 10)
fig1B
invisible(dev.off())


# -----------------------------------------------------------------------
# RNA v ATAC CORRELATIONS
# -----------------------------------------------------------------------

# common genes
common_genes <- intersect(Features(mean_expression), Features(mean_geneactivity))

# get scaled data matrix and order
df_rna <- LayerData(mean_expression, layer = "scale.data", features = common_genes)[common_genes, celltype_order$abbrev]
df_atac <- LayerData(mean_geneactivity, layer = "scale.data", features = common_genes)[common_genes, celltype_order$abbrev]
identical(rownames(df_rna), rownames(df_atac))

# plot correlations
df_correlation <- cor(df_rna, df_atac, method = "pearson")
identical(rownames(df_correlation), colnames(df_correlation))
identical(rownames(df_correlation), celltype_order$abbrev)

rownames(df_correlation) <- celltype_order$celltype
colnames(df_correlation) <- celltype_order$celltype

col_heatmap <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

fig1D <- Heatmap(df_correlation, cluster_rows = F, cluster_columns = F, 
                width = ncol(df_correlation)*unit(1,"cm"), 
                height = nrow(df_correlation)*unit(1,"cm"),
                row_names_side = "left", column_title_side = "bottom", col = col_heatmap,
                row_title = "RNA",
                column_title = "Gene activity",
                row_title_gp = gpar(fontsize = 18),
                column_title_gp = gpar(fontsize = 18),
                column_names_rot = 45,
                heatmap_legend_param = list(title = NULL, at = c(-1,0,1), labels = c(-1,0,1)))
fig1D

pdf(file.path(dir_out, "fig1D.pdf"), width = 6, height = 6)
fig1D
invisible(dev.off())


# correlation statistics for matched vs non-matched celltypes
diag_vals <- diag(df_correlation)
mask <- row(df_correlation) != col(df_correlation)
off_vals <- df_correlation[mask]
# and excl inh v exc neu comparisons
mask[7,8] <- FALSE
mask[8,7] <- FALSE
off_vals_excl_neu <- df_correlation[mask]

diag_stats <- calc_summary_stats(diag_vals)
off_stats <- calc_summary_stats(off_vals)
off_excl_neu_stats <- calc_summary_stats(off_vals_excl_neu)

write.table(data.frame(diag_entries = diag_stats, off_entries = off_stats, off_excl_neu_entries = off_excl_neu_stats),
            file = file.path(dir_out, "atac_rna_correlation_stats.tsv"), sep = "\t")

