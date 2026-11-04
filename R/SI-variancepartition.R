# Supplementary Figures/Tables for Figure 2 (variance partition results)

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PARAMS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(argparse)
  library(ComplexHeatmap)
  library(pals)
})

# parse arguments
parser <- ArgumentParser(description = "snMultiome figures")
parser$add_argument(
  "--config", required = TRUE, type = "character")

args <- parser$parse_args()

# config for file paths
source(args$config)


# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

dir_in <- file.path(dir_data, "/15_manuscript_figures/data")
dir_out <- file.path(dir_scripts, "data/figures-supplementary/variance-partitioning")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# variance partition results
vp_stack <- readRDS(file.path(dir_manuscript, "vp_dreamlet_stacked.rds"))
vp_celltype <- readRDS(file.path(dir_manuscript, "vp_dreamlet_celltype.rds"))

# cell type marker genes
marker_genes <- c("FGFR3", "AQP4", "P2RY12", "PLP1", "VCAN", "ZEB1", "SNAP25", "SYT1", "GAD1", "XIST", "UTY")

# named colours (for cell type ordering)
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# vp CCA results
vp_cca <- readRDS(file.path(dir_manuscript, "vp_dreamlet_cca.rds"))

# -----------------------------------------------------------------------
# VARIANCE PARTITION RESULTS - STACKED ASSAY
# -----------------------------------------------------------------------

# check total number of genes
sprintf("Total number of genes analysed: %s", length(unique(vp_stack$gene)))
# [1] "Total number of genes analysed: 11666"

# update variable names
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

# Variance explained by covariate
vp_stacked_stats <- vp_stack |>
  group_by(covariate) |>
  summarise(mean = mean(100 * variance),
            median = median(100 * variance),
            sd = sd(100 * variance)) |>
  arrange(desc(mean))

covar_order <- vp_stacked_stats |>
  dplyr::filter(covariate != "Residuals") |>
  pull(covariate)

vp_stacked_stats <- vp_stacked_stats |>
  mutate(covariate = factor(covariate, levels = c(covar_order, "Residuals"))) |>
  arrange(covariate)

write.csv(vp_stacked_stats, file.path(dir_out, "table_vp_stacked.csv"),
          row.names = F)


# -----------------------------------------------------------------------
# FIGURE SX: VP PER GENE
# -----------------------------------------------------------------------

# VP % for selected genes
vp_genes <- vp_stack |>
  filter(gene %in% marker_genes) |>
  mutate(gene = factor(gene, levels = marker_genes),
         Covariate = factor(covariate, levels = c(covar_order, "Residuals")))

# set colours for categories
cols_covariates <- setNames(c(cols25(length(covar_order)), "grey70"), 
                            nm = c(covar_order, "Residuals"))

# plot proportions
plot_vp <- vp_genes |>
  ggplot(aes(x = variance, y = gene, fill = Covariate)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = cols_covariates) +
  theme_minimal() +
  scale_y_discrete(limits = rev(marker_genes)) +
  xlab("Total variance explained (%)") + ylab("") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_line(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line.x = element_line(),
        panel.grid = element_blank(),
        legend.key.spacing.y = unit(2, "pt"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100))
  
ggsave(file.path(dir_out, "fig_vp_by_gene.pdf"), plot_vp, height = 4, width = 7)



# -----------------------------------------------------------------------
# VARIANCE PARTITION RESULTS - CELLTYPE ASSAY
# -----------------------------------------------------------------------

vp_celltype <- vp_celltype |>
  mutate(celltype = factor(assay)) 

# update cell type names
levels(vp_celltype$celltype) <- data.frame(celltype = levels(vp_celltype$celltype)) |>
  mutate(celltype = tools::toTitleCase(gsub("-", " ", celltype))) |>
  mutate(celltype = recode(celltype, "Oligodendrocyte Precursor" = "OPC")) |>
  pull(celltype)
vp_celltype <- vp_celltype |>
  mutate(celltype = factor(celltype, levels = names(named_colours$celltype_name))) |>
  select(!assay)

# update covariate names
vp_celltype <- vp_celltype |>
  tidyr::pivot_longer(cols = !c(gene, celltype), names_to = "covariate", values_to = "variance") |>
  mutate(covariate = replace_values(stringr::str_to_title(covariate),
                                    "Celltype" ~ "Cell type",
                                    "Scale.pmd" ~ "PMD (scaled)",
                                    "Scale.age" ~ "Age (scaled)",
                                    "Brainbank" ~ "Brain bank",
                                    "Log.nfeature_rna" ~ "No. of genes (log)",
                                    "Percent.mt" ~ "% mitochondrial",
                                    "Percent.ribosomal" ~ "% ribosomal"))

# Mean variance explained by covariate (per cell type)
vp_celltype_stats <- vp_celltype |>
  summarise(mean = mean(100 * variance),
            median = median(100 * variance),
            sd = sd(100 * variance), .by = c(celltype, covariate)) |>
  mutate(covariate = factor(covariate, levels = c(covar_order[-1], "Residuals"))) |>
  arrange(covariate, celltype)

write.csv(vp_celltype_stats, file.path(dir_out, "table_vp_celltype.csv"),
          row.names = F)


# -----------------------------------------------------------------------
# FIGURE SX: VP COVARIATE CORRELATIONS
# -----------------------------------------------------------------------

# CCA heatmap
covar_names <- data.frame(orig = rownames(vp_cca)) |>
  mutate(covariate = replace_values(stringr::str_to_title(orig),
                                    "Celltype" ~ "Cell type",
                                    "Scale.pmd" ~ "PMD (scaled)",
                                    "Scale.age" ~ "Age (scaled)",
                                    "Brainbank" ~ "Brain bank",
                                    "Log.nfeature_rna" ~ "No. of genes (log)",
                                    "Percent.mt" ~ "% mitochondrial",
                                    "Percent.ribosomal" ~ "% ribosomal")) |>
  mutate(covariate = factor(covariate), levels = covar_order)
covar_names <- setNames(covar_names$covariate, nm = covar_names$orig)

rownames(vp_cca) <- covar_names[rownames(vp_cca)]
colnames(vp_cca) <- covar_names[colnames(vp_cca)]

col_fun <- circlize::colorRamp2(c(0,1), c("white", "blue"))
p <- Heatmap(vp_cca,
             col = col_fun,
             column_names_rot = 45,
             width = nrow(vp_cca) * unit(1, "cm"),
             height = nrow(vp_cca) * unit(1, "cm"),
             heatmap_legend_param = list(title = "Correlation",
                                         title_gp = gpar(fontsize = 12),
                                         labels_gp = gpar(fontsize = 12)),
             column_names_gp = gpar(fontsize = 14),
             row_names_gp = gpar(fontsize = 14))
pdf(file.path(dir_out, "/fig_heatmap_cca.pdf"), width = 8, height = 7)
draw(p)
invisible(dev.off())

