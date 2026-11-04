# Functions for analysis-signac/R scripts

# figure-4.R

# -----------------------------------------------------------------------
# COVERAGE PLOTS
# -----------------------------------------------------------------------

plot_coverage_region <- function(sc, region, gene = gene, snp_gr, padding = 5000, resizing = 51){
  
  p <- CoveragePlot(
    object = sc,
    region = region,
    assay = "AMPFRAG",
    expression.assay = "RNA",
    features = gene,
    group.by = "celltype_name",
    region.highlight = resize(snp_gr, width = resizing, fix = "center"),
    extend.upstream = padding,
    extend.downstream = padding,
    heights = c(6, 1),
    peaks = FALSE,
    links = FALSE
  )
  
  return(p)
  
}

# option with more space on Genes track (for multiple genes in SI figures)
plot_coverage_region_SI <- function(sc, region, gene = gene, snp_gr, padding = 5000, resizing = 51){
  
  p <- CoveragePlot(
    object = sc,
    region = region,
    assay = "AMPFRAG",
    expression.assay = "RNA",
    features = gene,
    group.by = "celltype_name",
    region.highlight = resize(snp_gr, width = resizing, fix = "center"),
    extend.upstream = padding,
    extend.downstream = padding,
    heights = c(3, 1),
    peaks = FALSE,
    links = FALSE
  )
  
  return(p)
  
}
