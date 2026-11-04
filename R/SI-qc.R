# SI plots for sample QC and composition

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
  library(ggalluvial)
  library(ggrepel)
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
dir_out <- file.path(dir_scripts, "data/figures-supplementary/sample-qc")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# sc data
sc <- readRDS(file.path(dir_manuscript, "sc_rna.rds"))
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# marker genes
marker_genes <- read.csv(file.path(dir_scripts, "data/metadata/full-marker-genes.csv"), 
                         stringsAsFactors = F)

# Reference taxonomy cluster names
ref_celltypes <- read.csv(file.path(dir_scripts, "data/metadata/siletti_clusters.csv"), 
                          header = F, stringsAsFactors = F)

# -----------------------------------------------------------------------
# DATA PREP
# -----------------------------------------------------------------------

# replace donor IDs with numbers, numbered in order of total cell count
donor_counts <- sc@meta.data |>
  summarise(.by = individual, n = n()) |>
  arrange(desc(n)) |>
  mutate(donorID = row_number())

sc$donorID <- factor(sc$individual, levels = donor_counts$individual)
levels(sc$donorID) <- donor_counts$donorID

# -----------------------------------------------------------------------
# FigSX: QC metrics
# -----------------------------------------------------------------------

qc_metrics <- data.frame(metric = c("nCount_RNA", "nFeature_RNA", 
                                    "percent.mt", "percent.ribosomal",
                                    "nCount_AMPFRAG", "nFeature_AMPFRAG",
                                    "nucleosome_signal", "TSS.enrichment"),
                         newnames = c("log(nCount RNA)", "log(nFeature RNA)",
                                      "% mitochondrial", "% ribosomal",
                                      "log(nCount ATAC)", "log(nFeature ATAC)",
                                      "Nucleosome signal", "TSS enrichment"))

qc_summary <- sc@meta.data |>
  select(any_of(c("donorID", qc_metrics$metric))) |>
  tidyr::pivot_longer(cols = any_of(qc_metrics$metric), 
                      names_to = "metric",
                      values_to = "value") |>
  left_join(qc_metrics, by = "metric") |>
  mutate(newnames = factor(newnames, levels = qc_metrics$newnames)) |>
  mutate(value_plot = ifelse(metric %in% c("nCount_RNA", "nFeature_RNA",
                                           "nCount_AMPFRAG", "nFeature_AMPFRAG"),
                             log10(value), value))


# plot QC metrics grouped by donor ID
fig_qc <- qc_summary |>
  ggplot(aes(x = donorID, y = value_plot)) +
  geom_violin(scale = "width",
              color = "grey70", 
              fill = "grey70") +
  facet_wrap(vars(newnames),
             scales = "free_y",
             ncol = 4,
             axes = "all",
             axis.labels = "all") +
  theme_classic() +
  xlab("Donor ID") + ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
# fig_qc
ggsave(file.path(dir_out, "fig_qc_metrics.pdf"), plot = fig_qc, width = 12, height = 6)


# -----------------------------------------------------------------------
# FigSX: BARPLOTS OF CELL TYPE NUMBERS / PROPORTIONS
# -----------------------------------------------------------------------

# set celltype_name levels to match total count per cell type
cell_counts <- sc@meta.data |>
  summarise(.by = celltype_name, n = n()) |>
  arrange(desc(n))

# set sampleID order by donorID then region
sample_order <- as.vector(t(outer(1:11, names(named_colours$region), paste))) 
sample_order <- gsub(" ", "_", sample_order)

df_annot <- sc@meta.data |>
  select(celltype_name, subcelltype, region, donorID) |>
  mutate(celltype_name = factor(celltype_name, levels = rev(cell_counts$celltype_name))) |>
  mutate(sampleID = paste0(donorID, "_", region)) |>
  mutate(sampleID = factor(sampleID, levels = sample_order))

  
# frequency of cell types per region (SI to Fig 1C)
fig_region_celltype <- plot_barplot_vertical(df_annot, named_colours, "region", "celltype_name", proportion = F) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex")) +
  scale_x_continuous(breaks = c(0, 20000, 40000), labels = c(0, 20, 40)) +
  xlab(bquote("Number of nuclei "(10^3))) + 
  ylab("Region")

fig_region_celltype_proportion <- plot_barplot_vertical(df_annot, named_colours, "region", "celltype_name", proportion = T) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)") + 
  ylab("Region")

# frequency of cell types per donor
fig_donor_celltype <- plot_barplot_vertical(df_annot, named_colours, "donorID", "celltype_name", proportion = F) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = c(0, 10000, 20000), labels = c(0, 10, 20), limits = c(0, 20000)) +
  xlab(bquote("Number of nuclei "(10^3))) + 
  ylab("Donor ID")

fig_donor_celltype_proportion <- plot_barplot_vertical(df_annot, named_colours, "donorID", "celltype_name", proportion = T) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)") + 
  ylab("Donor ID")

# frequency of cell types per donor
fig_sample_celltype <- plot_barplot_vertical(df_annot, named_colours, "sampleID", "celltype_name", proportion = F) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = c(0, 4000, 8000), labels = c(0, 4, 8)) +
  xlab(bquote("Number of nuclei "(10^3))) + 
  ylab("Sample ID (donor / region)")

fig_sample_celltype_proportion <- plot_barplot_vertical(df_annot, named_colours, "sampleID", "celltype_name", proportion = T) +
  scale_fill_manual(NULL, values = named_colours$celltype_name,
                    breaks = names(named_colours$celltype_name)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)") + 
  ylab("Sample ID (donor / region)")

# Combined plot of sample cell type composition
fig_region_celltype + plot_spacer() + fig_region_celltype_proportion + 
  fig_donor_celltype + plot_spacer() + fig_donor_celltype_proportion +
  fig_sample_celltype + plot_spacer() + fig_sample_celltype_proportion +
  plot_layout(guides = "collect", axes = "collect", ncol = 3,
              heights = c(1, 11/4, 11), 
              widths = c(1.5, 0.1, 1))
ggsave(file.path(dir_out, "figSX_celltype_composition.pdf"), width = 10, height = 12)



# -----------------------------------------------------------------------
# FigSX: BARPLOTS OF SUB CELL TYPE NUMBERS / PROPORTIONS
# -----------------------------------------------------------------------

# region x subcelltype
fig_region_subcelltype <- plot_barplot_vertical(df_annot, named_colours, "region", "subcelltype", proportion = F) +
  scale_fill_manual(NULL, values = named_colours$subcelltype, breaks = names(named_colours$subcelltype)) +
  scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex"))  +
  scale_x_continuous(breaks = c(0, 20000, 40000), labels = c(0, 20, 40)) +
  xlab(bquote("Number of nuclei "(10^3))) + 
  ylab("Region") +
  NoLegend()

fig_region_subcelltype_proportion <- plot_barplot_vertical(df_annot, named_colours, "region", "subcelltype", proportion = T) +
  scale_fill_manual(NULL, values = named_colours$subcelltype, breaks = names(named_colours$subcelltype)) +
  scale_y_discrete(limits = rev, labels = c("PFC" = "Prefrontal cortex"))  +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)") + 
  ylab("Region") +
  theme(legend.key.height = unit(12, "pt"),
        legend.key.width = unit(12, "pt"),
        legend.key.spacing.y = unit(1, "pt"),
        legend.text = element_text(size = 11))

fig_region_subcelltype + plot_spacer() + fig_region_subcelltype_proportion + 
  plot_layout(guides = "collect", axes = "collect", widths = c(1, 0.1, 1))
ggsave(file.path(dir_out, "figSX_region_subcelltype_composition.pdf"), width = 12, height = 3)

# subcelltype x region (only include neurons or axes skewed by high numbers of oligodendrocytes)
named_colours$Region <- setNames(named_colours$region,
                                 nm = gsub("PFC", "Prefrontal cortex", names(named_colours$region)))
df_neuron <- sc@meta.data |>
  select(subcelltype, region, cellclass) |>
  filter(cellclass == "neuronal") |>
  mutate(Region = recode(region, "PFC" = "Prefrontal cortex")) |>
  droplevels()
  
neuron_colours <- list(subcelltype = named_colours$subcelltype[-c(1:6)], Region = named_colours$Region)
fig_neurons_region <- plot_barplot_vertical(df_neuron, neuron_colours,
                                          "subcelltype", "Region", proportion = F) + 
  scale_x_continuous(breaks = c(0, 5000, 10000), labels = c(0, 5, 10)) +
  xlab(bquote("Number of nuclei "(10^3))) + 
  ylab("Neuron subtype") +
  NoLegend()

fig_neurons_region_proportion <- plot_barplot_vertical(df_neuron, neuron_colours,
                                            "subcelltype", "Region") + 
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  xlab("Proportion of nuclei (%)")

fig_neurons_region + plot_spacer() + fig_neurons_region_proportion + 
  plot_layout(guides = "collect", axes = "collect", widths = c(1, 0.1, 0.8))
ggsave(file.path(dir_out, "figSX_neuron_region_composition.pdf"), width = 10, height = 6)

# record stats on regional specificity
df_neu_stats <- df_neuron |>
  count(subcelltype, Region, name = "n") |>
  mutate(n_subtype = sum(n),
         pc_subtype = n / n_subtype,
         max_pc = max(pc_subtype),
         .by = subcelltype) |>
  arrange(desc(max_pc), subcelltype)

write.csv(df_neu_stats, file = file.path(dir_out, "/table_neu_subtype_regions.csv"),
          row.names = F)


# -----------------------------------------------------------------------
# FigSX: COMPARISON TO REFERENCE ANNOTATIONS - ALLUVIAL
# -----------------------------------------------------------------------

# Alluvial map of cell types vs reference annotations
df_compare <- sc@meta.data |>
  group_by(celltype_name, allen_supercluster) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(freq = n/sum(n)) |>
  mutate(allen_supercluster = factor(allen_supercluster, levels = ref_celltypes$V1),
         celltype_name = factor(tolower(celltype_name), levels = tolower(names(named_colours$celltype_name))))

# converted celltype_name to lower case to avoid overlapping names vs reference taxonomy
# use helper function to format correctly for text labels in plot
format_text_label <- function(x) {
  x <- gsub("opc", "OPC", x)
  x <- tools::toTitleCase(x)
  return(x)
}

alluvial_colours <- setNames(named_colours$celltype_name, nm = tolower(names(named_colours$celltype_name)))
                             
fig_alluvial <- df_compare |>
  ggplot(aes(axis1 = allen_supercluster, axis2 = celltype_name, y = freq)) +
  geom_alluvium(aes(fill = celltype_name)) +
  geom_stratum(color = "grey70", aes(alpha = 0.5)) +
  # label cell types in stratum
  geom_text(stat = "stratum", aes(label = ifelse(after_stat(x) == 2, format_text_label(as.character(after_stat(stratum))), "")), size = 5) +
  geom_text(stat = "stratum", aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")), size = 5, min.y = 0.02) +
  # repel labels where needed for c1
  geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")), size = 5, nudge_x = -0.7, min.y = 0.005, max.y = 0.02) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.15))) +
  theme_minimal() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  # add labels
  annotate("text", x = c(1,2), y = c(1.04, 1.04), label = c("Reference taxonomy", "Annotation"), size = 6) +
  labs(x = NULL) +
  scale_fill_manual(values = alluvial_colours) +
  ylab("Proportion of nuclei (%)") +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    axis.ticks.y = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100))

fig_alluvial
ggsave(file.path(dir_out, "figSX_alluvial.pdf"), fig_alluvial, width = 10, height = 8)


# -----------------------------------------------------------------------
# FigSX: COMPARISON TO REFERENCE ANNOTATIONS - HEATMAP
# -----------------------------------------------------------------------

# % of cells within each of manual subcell type annotations that fall within allen superclusters
df_subtypes <- sc@meta.data |>
  mutate(allen_supercluster = factor(allen_supercluster, levels = ref_celltypes$V1)) |>
  group_by(subcelltype, allen_supercluster) |>
  summarise(n = n(), .groups = "drop")
mat_subtypes <- xtabs(n ~ allen_supercluster + subcelltype, data = df_subtypes) 
mat_subtypes <- proportions(mat_subtypes, margin = 2) * 100

fig_heatmap <- Heatmap(mat_subtypes, cluster_rows = F, cluster_columns = F,
                       col = c("white", "blue"),
                       width = dim(mat_subtypes)[2] * unit(0.8, "cm"),
                       height = dim(mat_subtypes)[1] * unit(0.8, "cm"),
                       row_names_side = "left",
                       column_names_rot = 45,
                       column_title_side = "bottom", column_title = "Cell type annotation",
                       column_title_gp = gpar(fontsize = 16),
                       column_names_gp = gpar(fontsize = 12),
                       row_title = "Reference taxonomy supercluster",
                       row_title_gp = gpar(fontsize = 16),
                       row_names_gp = gpar(fontsize = 12),
                       row_names_max_width = unit(7,"cm"),
                       heatmap_legend_param = list(
                         title = "% nuclei",
                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                         labels_gp = gpar(fontsize = 12)),
                       cell_fun = function(j, i, x, y, width, height, fill){
                         cellval <- mat_subtypes[i,j]
                         if (cellval > 1 & cellval < 50) {
                           grid.text(sprintf("%.0f", cellval), x, y, gp = gpar(fontsize = 10))
                         } else if (cellval > 50) {
                           grid.text(sprintf("%.0f", cellval), x, y, gp = gpar(fontsize = 10, col = "white"))
                         }})

pdf(file.path(dir_out, "figSX_annotation_comparison_heatmap.pdf"), width = 12, height = 12)
fig_heatmap
invisible(dev.off())


