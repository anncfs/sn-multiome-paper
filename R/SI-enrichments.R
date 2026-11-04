# SI figures (heritability / enrichments)

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
dir_out <- file.path(dir_scripts, "data/figures-supplementary/enrichments")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# EWCE results
df_ewce <- readRDS(file.path(dir_manuscript, "ewce_results_level1_10000.rds"))

# colours
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# combined heritability enrichment results
df_heritability <- readRDS(file.path(dir_manuscript, "heritability_enrichment_results.rds")) |>
  mutate(CELLTYPE = gsub("-", "_", CELLTYPE))

# GWAS info
df_gwas <- read.csv(file.path(dir_manuscript, "gwas_info.csv"), header = T) |>
  as_tibble() |>
  mutate(Category = factor(Category, levels = names(named_colours$Category))) |>
  arrange(Category, GWAS)
df_gwas$GWAS <- factor(df_gwas$GWAS, levels = df_gwas$GWAS)

# gene specificities
df_specificity <- read.csv(file.path(dir_manuscript, "gene_properties_subcelltype.txt"), sep = "\t", stringsAsFactors = F)

# -----------------------------------------------------------------------
# PLOT EWCE RESULTS
# -----------------------------------------------------------------------

# update cell type labels
df_ewce <- df_ewce |>
  mutate(celltype = tools::toTitleCase(gsub("_", " ", CellType))) |>
  mutate(celltype = recode(celltype, "Oligodendrocyte Precursor" = "OPC")) |>
  mutate(celltype = factor(celltype, levels = names(named_colours$celltype_name)),
         gene_set = factor(gene_set)) |>
  mutate(gene_window = recode(gene_window,
                              "10_2_10000" = "-10/+2 kb",
                              "35_10_10000" = "-35/+10 kb"))

maxeffect <- ceiling(max(abs(df_ewce$sd_from_mean)))

p_list <- lapply(unique(df_ewce$gene_window), 
                 function(x) plot_ewce_enrichment(filter(df_ewce, gene_window == x),
                                                  maxeffect = maxeffect))

lgd_list <- list(Legend(title = "FDR", labels = c("< 0.05", "< 0.01"),
       title_gp = gpar(fontsize = 12, fontface = "bold"),
       labels_gp = gpar(fontsize = 12),
       graphics = list(function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 14)),
                       function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 14)))))

pdf(file.path(dir_out, "ewce_results.pdf"), width = 9, height = 5)
draw(p_list[[1]] + p_list[[3]], ht_gap = unit(1, "cm"), merge_legends = TRUE, annotation_legend_list = lgd_list)
dev.off()


# -----------------------------------------------------------------------
# Fig SX: MAGMA MAJOR CELLTYPE
# -----------------------------------------------------------------------

# unconditioned results, major celltype, all GWAS
df <- df_heritability |>
  filter(ANALYSIS == "unconditioned", GROUPING == "celltype", GWAS %in% df_gwas$GWAS) |>
  mutate(CELLTYPE = factor(CELLTYPE))

# udpate celltype names and order
levels(df$CELLTYPE) <- sapply(seq_along(levels(df$CELLTYPE)), function(x) {
  grep(levels(df$CELLTYPE)[x], names(named_colours$celltype_name), value = T)
})
df <- df |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$celltype_name)))

pdf(file.path(dir_out, "fig_magma_celltype.pdf"), width = 15, height = 6)
plot_heatmap_enrichment(df, df_gwas = df_gwas)
dev.off()


# -----------------------------------------------------------------------
# Fig SX: BMI CONDITIONAL ANALYSIS
# -----------------------------------------------------------------------

df <- df_heritability |>
  filter(METHOD == "MAGMA", GWAS == "BMI_2018", GROUPING == "subcelltype") |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$subcelltype))) |>
  mutate(CONDITIONED_ON = gsub("-", "_", CONDITIONED_ON))

max_beta <- ceiling(max(abs(df$BETA)))

# iterative plot
p_iter <- plot_heatmap_conditional_iterative(df |> filter(ANALYSIS == "iterative"), 
                                             id = "BMI_2018", title = "Iterative")
p_iter@heatmap_param$show_heatmap_legend <- FALSE

# per cell type plot
p_percelltype <- plot_heatmap_conditional_percelltype(df |> filter(ANALYSIS == "each"),
                                                      id = "BMI_2018",
                                                      title = "Per celltype",
                                                      named_colours = named_colours)

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

pdf(file.path(dir_out, "fig_bmi_conditional_magma.pdf"), width = 12, height = 10)
draw(p_iter + p_percelltype, ht_gap = unit(1, "cm"), merge_legend = TRUE, annotation_legend_list = lgd_list)
dev.off()


# -----------------------------------------------------------------------
# Fig SX: GENE SPECIFICITY CORRELATIONS
# -----------------------------------------------------------------------

colnames(df_specificity) <- gsub("\\.", "_", names(df_specificity))
corr_specificity <- df_specificity |>
  select(names(named_colours$subcelltype)) |>
  cor(method = "pearson")

fig_corr <- Heatmap(corr_specificity,
                    width = dim(corr_specificity)[2] * unit(0.8, "cm"),
                    height = dim(corr_specificity)[1] * unit(0.8, "cm"),
                    column_names_rot = 45,
                    column_names_gp = gpar(fontsize = 14),
                    row_names_gp = gpar(fontsize = 14),
                    heatmap_legend_param = list(title = "Correlation",
                                                at = seq.int(-1, 1, 1),
                                                labels_gp = gpar(fontsize = 14),
                                                title = gpar(fontsize = 14, fontface = "bold")))

pdf(file.path(dir_out, "fig_specificity_correlations.pdf"), width = 13, height = 10)
fig_corr
dev.off()



# -----------------------------------------------------------------------
# Fig SX: LDSC RESULTS
# -----------------------------------------------------------------------

# restrict to main figure GWAS
df_gwas <- df_gwas |>
  filter(Plot == TRUE) |>
  arrange(Category, GWAS)
df_gwas$GWAS <- factor(df_gwas$GWAS, levels = df_gwas$GWAS)

# LDSC RNA analysis
df <- df_heritability |>
  filter(METHOD == "LDSC", GROUPING == "celltype") |>
  mutate(CELLTYPE = factor(CELLTYPE)) |>
  filter(GWAS %in% df_gwas$GWAS,
         MODALITY == "RNA")

# udpate celltype names and order
levels(df$CELLTYPE) <- sapply(seq_along(levels(df$CELLTYPE)), function(x) {
  grep(levels(df$CELLTYPE)[x], names(named_colours$celltype_name), value = T)
})
df <- df |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$celltype_name)))

fig_ldsc_rna <- plot_heatmap_enrichment(df = df,
                                        coefficient_name = "TAU.NORM", 
                                        df_gwas = df_gwas,
                                        colour = "darkgreen")


# LDSC ATAC SUBCELL TYPE
df <- df_heritability |>
  filter(GROUPING == "subcelltype", MODALITY == "ATAC", GWAS %in% df_gwas$GWAS) |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$subcelltype)))

fig_ldsc_subcelltype <- plot_heatmap_enrichment(df = df,
                                                coefficient_name = "TAU.NORM",
                                                df_gwas = df_gwas,
                                                colour = "royalblue4")

pdf(file.path(dir_out, "figSX_LDSC_results.pdf"), width = 13, height = 14)
fig_ldsc_rna@bottom_annotation <- NULL
draw(fig_ldsc_rna %v% fig_ldsc_subcelltype,
     merge_legend = FALSE, ht_gap = unit(1.5, "cm"),
     legend_grouping = "original")
dev.off()



# -----------------------------------------------------------------------
# SI TABLES
# -----------------------------------------------------------------------

# reload GWAS info as filtered above
df_gwas <- read.csv(file.path(dir_manuscript, "gwas_info.csv"), header = T) |>
  as_tibble() |>
  mutate(Category = factor(Category, levels = names(named_colours$Category))) |>
  arrange(Category, GWAS)
df_gwas$GWAS <- factor(df_gwas$GWAS, levels = df_gwas$GWAS)

# Major cell type results, unconditioned --------------------------
df_celltype <- df_heritability |>
  filter(GROUPING == "celltype", GWAS %in% df_gwas$GWAS, MODALITY %in% c("ATAC", "RNA")) |>
  mutate(CELLTYPE = factor(CELLTYPE))

# udpate celltype names and order
levels(df_celltype$CELLTYPE) <- sapply(seq_along(levels(df_celltype$CELLTYPE)), function(x) {
  grep(levels(df_celltype$CELLTYPE)[x], names(named_colours$celltype_name), value = T)
})
df_celltype <- df_celltype |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$celltype_name))) |>
  arrange(desc(METHOD), desc(MODALITY), GWAS, CELLTYPE) |>
  select(MODALITY, GWAS, CELLTYPE, BETA, TAU.NORM, FDR, METHOD)


# save MAGMA & sLDSC cell type results
write.table(df_celltype, file = file.path(dir_out, "tableSX_celltype_enrichments.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)




# Sub cell type results, unconditioned --------------------------
df_subcelltype <- df_heritability |>
  filter(GROUPING == "subcelltype", GWAS %in% df_gwas$GWAS, 
         MODALITY %in% c("ATAC", "RNA"), !(ANALYSIS %in% c("iterative", "each"))) |>
  mutate(CELLTYPE = factor(CELLTYPE, levels = names(named_colours$subcelltype))) |>
  arrange(desc(METHOD), desc(MODALITY), GWAS, CELLTYPE) |>
  select(MODALITY, GWAS, CELLTYPE, BETA, TAU.NORM, FDR, METHOD)

write.table(df_subcelltype, file = file.path(dir_out, "tableSX_subcelltype_enrichments.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
