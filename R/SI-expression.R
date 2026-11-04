# SI plots for marker gene expression

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
dir_out <- file.path(dir_scripts, "data/figures-supplementary/gene-expression")
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

# -----------------------------------------------------------------------
# FigSX: MARKER GENE EXPRESSION
# -----------------------------------------------------------------------

# mean gene expression (neurons)
gene_order <- marker_genes |>
  filter(grepl("Neuron", Celltype)) |>
  pull(Gene)

# neurons only
celltype_order <- setdiff(levels(sc$subcelltype), levels(sc$celltype_name))

# Idents(sc) <- "subcelltype"
DotPlot(sc, features = gene_order, group.by = "subcelltype", 
         cols = c("grey90", "blue"),scale = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 16)) +
  ylim(celltype_order) +
  ylab("Neuron subtype") + xlab("Gene")
ggsave(file.path(dir_out, "fig_neuron_expression.pdf"), width = 10, height = 6)
  
