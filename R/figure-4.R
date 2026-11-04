# Figure 4 plots (peak-gene links)

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PARAMS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
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

# functions
source(file.path(dir_scripts, "analysis/analysis-signac/R/functions-signac.R"))


# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

dir_in <- paste0(dir_data, "/15_manuscript_figures/data")
dir_out <- file.path(dir_scripts, "data/figures-main/figure-4")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

dir_SI <- file.path(dir_scripts, "data/figures-supplementary/snp-gene-links")
if(!file.exists(dir_SI)){dir.create(dir_SI, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# load multiome data
sc_rna <- readRDS(paste0(dir_in, "/sc_rna.rds"))
sc_atac <- readRDS(paste0(dir_in, "/sc_atac.rds"))
sc_rna[["AMPFRAG"]] <- sc_atac[["AMPFRAG"]]

# colours
named_colours <- readRDS(file.path(dir_manuscript, "named_colour_schemes.rds"))

# AD credible set lead SNPs
df_ad <- readRDS(paste0(dir_in, "/df_pgboost_links_ad_credible_set.rds"))

# pgBoost links
df_links <- readRDS(paste0(dir_in, "/df_pgboost_links.rds"))

# -----------------------------------------------------------------------
# FORMAT DATA
# -----------------------------------------------------------------------

# add region stats info for coverage plotting
DefaultAssay(sc_rna) = "AMPFRAG"
sc_rna <- RegionStats(sc_rna, genome = BSgenome.Hsapiens.UCSC.hg38)

# -----------------------------------------------------------------------
# PLOT COVERAGE PLOTS - INDIVIDUAL
# -----------------------------------------------------------------------

DefaultAssay(sc_rna) = "RNA"
Idents(sc_rna) <- sc_rna$celltype_name

# # re order cells
# named_colours$celltype_name <- named_colours$celltype_name[c(1:2,5:8,3:4)]
# sc_rna$celltype_name <- factor(sc_rna$celltype_name, levels = named_colours$celltype_name)

# selected SNPs and genes
list_snps <- list(
  RASGEF1C = c(rs113706587 = "chr5:180201150"),
  PICALM = c(rs10792832 = "chr11:86156833", rs10792831 = "chr11:86141937"),
  INPP5D = c(rs10933431 = "chr2:233117202")
)

# plot coverage for each gene + linked snps
p_list <- list()
for (gene in names(list_snps)) {

  # granges for snps
  snps_gr <- GRanges(list_snps[[gene]], 
    color = "grey10", 
    name = names(list_snps[[gene]])
  )

  # center and width for plot
  midpoint <- GRanges(seqnames(snps_gr[1]), mean(start(snps_gr)))
  width <- max(5000, start(midpoint) - min(start(snps_gr)))
  width <- ceiling(width/1000) * 1000

  p_list[[gene]] <- plot_coverage_region(sc_rna, region = midpoint, 
    gene = gene, snp_gr = snps_gr, padding = width)

  pdf(file.path(dir_out, paste0("/coverage_plot_", gene, ".pdf")), width = 6, height = 4)
  print(p_list[[gene]] & 
    ggplot2::scale_fill_manual(values = named_colours$celltype_name) &
    ggplot2::theme(text = ggplot2::element_text(size = 12),
      strip.text.x = ggplot2::element_text(hjust = 0.5)))
  invisible(dev.off())
}


# -----------------------------------------------------------------------
# PLOT COVERAGE PLOTS - ALL TOP AD CREDIBLE SET SNPs (SI)
# -----------------------------------------------------------------------

# format list of SNPs / genes (all AD lead SNPs with pgboost brain score > 0.95)
snps_ad <- df_ad |>
  select(snp, gene) |>
  distinct() |>
  tidyr::separate_wider_delim(cols = snp, delim = "-", names = c("chr", "pos-1", "bp")) |>
  select(!'pos-1') |>
  mutate(snp = paste0(chr, ":", bp))

# plot each linked SNP-gene pair
pdf(file.path(dir_SI, "/coverage_plots_ad_lead_snps.pdf"), width = 6, height = 4)

for (i in 1:nrow(snps_ad)) {
# for (i in 3) {
  gene <- snps_ad$gene[i]
  snp <- snps_ad$snp[i]

  # granges for snps
  snp_gr <- GRanges(snp, 
    color = "grey10", 
    name = snp
  )

  # plot centered on snp
  p <- plot_coverage_region_SI(sc_rna, region = snp_gr, 
    gene = gene, snp_gr = snp_gr, padding = 5000)

  print(p & 
    ggplot2::scale_fill_manual(values = named_colours$celltype_name) &
    ggplot2::theme(text = ggplot2::element_text(size = 12)))
  
}

invisible(dev.off())

# table of AD credible set lead SNP details
write.table(df_ad, file.path(dir_SI, "ad_credible_set_lead_snps.tsv"), sep = "\t", row.names = F, quote = F)

# average expression of AD linked genes
mean_ad_genes <- AverageExpression(sc_rna, assays = "RNA", features = unique(df_ad$gene), group.by = "celltype_name")
write.table(mean_ad_genes$RNA, file.path(dir_SI, "ad_credible_set_linked_genes_average_expression.tsv"), sep = "\t", row.names = T, quote = F)

