# Functions for manuscript figure plotting scripts

# figure-1.R
# figure-2.R
# figure-3.R
# SI-enrichments.R
# SI-qc.R

# -----------------------------------------------------------------------
# GENE EXPRESSION HEATMAP
# -----------------------------------------------------------------------

# heatmap of mean gene expression / gene activity
get_heatmap_expression <- function(dfX, col_limit = 2.5, title = "", legend_title = "") {
  col_heatmap <- circlize::colorRamp2(c(-col_limit, 0, col_limit), c("blue", "white", "red"))
  pX <- Heatmap(dfX, cluster_rows = F, cluster_columns = F, 
                width = ncol(dfX)*unit(0.7,"cm"), height = nrow(dfX)*unit(1,"cm"),
                row_names_side = "left", col = col_heatmap,
                # column_title = title,
                # column_title_gp = gpar(fontsize = 20),
                row_names_gp = gpar(fontsize = 14),
                column_names_rot = 45,
                heatmap_legend_param = list(title = legend_title, 
                                            # direction = "horizontal",
                                            title_position = "topcenter",
                                            at = c(-col_limit,0,col_limit), labels = c(-col_limit, 0, col_limit)),
                top_annotation = HeatmapAnnotation(
                  lab = anno_block(labels = title,
                                   gp = gpar(fill = NA, col = NA),
                                   labels_gp = gpar(fontsize = 20)),
                  annotation_height = unit(8, "mm"),
                  show_annotation_name = FALSE
                ))
  return(pX)
}



# -----------------------------------------------------------------------
# VERTICAL BAR PLOTS (FREQ / PROPORTION)
# -----------------------------------------------------------------------

# proportion or frequency of cells per category, vertical
plot_barplot_vertical <- function(df_annot, named_colours, id_group, id_colour,
                                              proportion = TRUE,
                                              group_colours = FALSE){
  
  # count cells per category
  tb <- df_annot |>
    group_by(.data[[id_group]], .data[[id_colour]]) |>
    summarise(ncells = n(), .groups = "drop") |>
    group_by(.data[[id_group]]) |>
    mutate(total = sum(ncells)) |>
    arrange(desc(total))
  
  # set y axis order
  y_order <- rev(names(named_colours[[id_group]]))
  
  # bar plot
  p1 <- tb %>%
    ggplot(aes(x = ncells, y = .data[[id_group]], fill = .data[[id_colour]]))
  
  # proportion v frequency
  if (proportion == TRUE) {
    p1 <- p1 + 
      geom_col(position = position_fill(), width = 0.8) +
      scale_x_continuous(breaks = seq(0, 1, 0.5), expand = expansion(mult = c(0, 0)))
  } else {
    p1 <- p1 +
      geom_col(width = 0.8)
  }

  p1 <- p1 +
    scale_y_discrete(limits = y_order) +
    theme_minimal() +
    scale_fill_manual(tools::toTitleCase(id_colour), values = named_colours[[id_colour]]) +
    ylab("") + xlab("") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.ticks.x = element_line(),
          axis.ticks.length = unit(0.2, "cm"),
          axis.line.x = element_line(),
          panel.grid = element_blank(),
          plot.margin = margin(4, 0, 4, 0),
          legend.key.spacing.y = unit(2, "pt"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold")) 
  
  # subcelltype colours
  p2 <- tb |>
    distinct(.data[[id_group]]) |>
    mutate(x = 1) |>
    ggplot(aes(x = x, y = .data[[id_group]], colour = .data[[id_group]])) +
      geom_point(size = 3) +
      scale_colour_manual(values = named_colours[[id_group]], guide = "none") +
      theme_minimal() +
      scale_y_discrete(limits = y_order) +
      scale_x_continuous(limits = c(0.7, 1.3), expand = expansion(mult = c(0, 0))) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.text.y.left = element_text(margin = margin(r = 6)),
            plot.margin = margin(4, 0, 4, 0))
  
  # plot layout
  if (group_colours) {
    p_combined <- p2 + p1 + theme(axis.text.y = element_blank()) + plot_layout(widths = c(0.08, 1))
  } else {
    p_combined <- p1
  }
 
  return(p_combined)
  
}


# -----------------------------------------------------------------------
# SUMMARY STATS (CORRELATION MATRIX)
# -----------------------------------------------------------------------

calc_summary_stats <- function(x) {
  c(min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE))
}


# -----------------------------------------------------------------------
# HERITABILITY HEATMAPS (UNCONDITIONED ANALYSES)
# -----------------------------------------------------------------------

plot_heatmap_enrichment <- function(df, maxeffect = NA, star = TRUE, coefficient_name = "BETA", 
                                    df_gwas = NULL, title = "", colour = "red"){

  # set colours
  if (is.na(maxeffect)) {
    maxeffect <- ceiling(max(abs(df[,coefficient_name])))
  }
  col_fun <- circlize::colorRamp2(c(-maxeffect, 0, maxeffect), c("grey", "white", colour))
  
  # plot ordering
  celltype_order <- levels(df$CELLTYPE)
  gwas_order <- levels(df_gwas$GWAS)
  
  # fdr matrix
  mat_fdr <- df |>
    select(GWAS, CELLTYPE, FDR) |>
    tidyr::pivot_wider(names_from = GWAS, values_from = FDR) |>
    tibble::column_to_rownames("CELLTYPE")
  mat_fdr <- mat_fdr[celltype_order, gwas_order] |>
    as.matrix()
  rownames(mat_fdr) <- celltype_order
  
  # effect size matrix
  mat_effect <- df %>%
    select(GWAS, CELLTYPE, all_of(coefficient_name)) %>%
    tidyr::pivot_wider(names_from = GWAS, values_from = all_of(coefficient_name)) %>%
    tibble::column_to_rownames("CELLTYPE")
  mat_effect <- mat_effect[celltype_order, gwas_order] |>
    as.matrix()
  rownames(mat_effect) <- celltype_order
  
  # get annotations
  plot_annotations <- add_heatmap_annotations(df_gwas, named_colours, celltype_order)
  
  
  coeff_label <- ifelse(coefficient_name == "TAU.NORM", "TAU*", coefficient_name)
  p <- Heatmap(mat_effect, 
               cluster_rows = F, cluster_columns = F, col = col_fun,
                 width = dim(mat_effect)[2] * unit(0.8, "cm"),
                 height = dim(mat_effect)[1] * unit(0.8, "cm"),
                 row_names_side = "left",
                 column_labels = df_gwas$Trait,
                 column_names_rot = 45,
                 column_title = title,
                 column_title_gp = gpar(fontsize = 16),
                 column_names_gp = gpar(fontsize = 16),
                 row_title_gp = gpar(fontsize = 16),
                 row_names_gp = gpar(fontsize = 16),
                 left_annotation = plot_annotations[[1]],
                 bottom_annotation = plot_annotations[[2]],
                 heatmap_legend_param = list(title = coeff_label, 
                                             title_gp = gpar(fontsize = 14, fontface = "bold"),
                                             at = c(-maxeffect, 0, maxeffect),
                                             labels_gp = gpar(fontsize = 14)),
                 border = T,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   testval <- mat_fdr[i, j]
                   if (star == TRUE & testval <= 0.01) {
                     grid.text("**", x, y, gp = gpar(fontsize = 14))
                   } else if (star == TRUE & testval <= 0.05) {
                     grid.text("*", x, y, gp = gpar(fontsize = 14))
                   }
                 })
  
  
  
  return(p)
}



# add heatmap annotations (cell types in rows, traits in columns)
add_heatmap_annotations <- function(df_gwas, named_colours, celltype_order) {
  
  # trait annotations
  ha_col <- columnAnnotation(df = select(df_gwas, Category),
                             col = named_colours,
                             annotation_legend_param = list(Category = list(title = "Category",
                                                                            title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                                            labels_gp = gpar(fontsize = 14))),
                             show_annotation_name = F,
                             show_legend = T)
  
  # cell type annotations
  col_neuron <- setNames(c("grey90", "grey60"), c(TRUE, FALSE))
  df_row <- data.frame(celltype = celltype_order) |>
    mutate(Neuron = ifelse(grepl("Exc|Inh|Neu", celltype), TRUE, FALSE))
  ha_row <- rowAnnotation(df = select(df_row, Neuron),
                          col = list(Neuron = col_neuron),
                          annotation_legend_param = list(Neuron = list(title = "Neuron",
                                                                       title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                                       labels_gp = gpar(fontsize = 14))),
                          show_annotation_name = F, show_legend = T)

  # # legends
  # lgd_col <- Legend(labels = names(named_colours$Category), title = "Trait category", type = "grid", legend_gp = gpar(fill = unname(named_colours$Category)))
  # lgd_row <- Legend(labels = names(col_neuron), title = "Neuron", type = "grid", legend_gp = gpar(fill = unname(col_neuron)))
  
  return(list(ha_row, ha_col))
}




# -----------------------------------------------------------------------
# MAGMA ENRICHMENTS (ITERATIVE ANALYSES)
# -----------------------------------------------------------------------

# if star = TRUE, then indicates significance of FDR (* = < 0.05, ** = < 0.01)

# heatmap of BETA/FDR for each trait, where row = celltype, and column indicates conditioned on celltype

plot_heatmap_conditional_iterative <- function(df, id = "BMI_2018", max_beta = NA, star = TRUE, title = NA,
                                               named_colours = named_colours) {
  
  # prep data
  dfx <- df %>%
    filter(GWAS == id) |>
    arrange(CELLTYPE, ROUND)
  
  # find last round where signif cell type remains
  last_signif_round <- dfx |>
    filter(FDR < 0.05) |>
    select(ROUND) |>
    max()
  
  # get order of conditioning labels (i.e. cell types) in each round
  # missing labels are those that were conditioned on in previous rounds, so not included in current round
  # new_conditioner <- dfx |>
  #   filter(ROUND <= last_signif_round + 1) |>
  #   select(ROUND, CELLTYPE) |>
  #   group_by(ROUND) |>
  #   summarise(tested = list(CELLTYPE), .groups = "drop") |>
  #   mutate(
  #     conditioned_on = lapply(tested, function(x) setdiff(levels(dfx$CELLTYPE), x)),
  #     conditioned_prev = lag(conditioned_on, default = list(character(0))),
  #     new_conditioner = mapply(setdiff, conditioned_on, conditioned_prev, SIMPLIFY = FALSE)
  #   ) |>
  #   select(ROUND, new_conditioner) |>
  #   tidyr::unnest(cols = new_conditioner)
  
  
  # set colours
  if (is.na(max_beta)) {
    max_beta <- ceiling(max(abs(dfx$BETA)))
  }
  col_fun <- circlize::colorRamp2(c(-max_beta, 0, max_beta), c("grey", "white", "red"))
  
  # create matrix (fdr and beta)
  mat_fdr <- dfx |>
    filter(ROUND <= last_signif_round + 1) |>
    select(CELLTYPE, ROUND, FDR) |>
    tidyr::pivot_wider(names_from = ROUND, values_from = FDR, values_fill = Inf)
  mat_fdr <- mat_fdr[,-1]
  
  mat_beta <- dfx |>
    filter(ROUND <= last_signif_round + 1) |>
    select(CELLTYPE, ROUND, BETA) |>
    tidyr::pivot_wider(names_from = ROUND, values_from = BETA, values_fill = list(value = NA_real_)) |>
    tibble::column_to_rownames(var = "CELLTYPE")  |>
    as.matrix()
  
  # create plot
  plt <- Heatmap(mat_beta, cluster_rows = F, cluster_columns = FALSE, col = col_fun,
                 width = dim(mat_beta)[2] * unit(0.7, "cm"),
                 height = dim(mat_beta)[1] * unit(0.8, "cm"),
                 row_names_side = "left",
                 na_col = "white",
                 column_title_gp = gpar(fontsize = 15),
                 column_names_gp = gpar(fontsize = 15),
                 row_names_gp = gpar(fontsize = 15),
                 column_title = ifelse(is.na(title), id, title),
                 # column_labels = c("-", new_conditioner$new_conditioner),
                 column_names_rot = 0,
                 column_names_centered = TRUE,
                 heatmap_legend_param = list(
                   title = "BETA",
                   at = seq.int(-max_beta, max_beta, max_beta),
                   labels = seq.int(-max_beta, max_beta, max_beta),
                   labels_gp = gpar(fontsize = 12),
                   title_gp = gpar(fontsize = 12, fontface = "bold")
                 ), border = T,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   # signifant results (*, **)
                   testval <- mat_fdr[i, j]
                   if (star == TRUE & testval <= 0.01) {
                     grid.text("**", x, y, gp = gpar(fontsize = 14))
                   } else if (star == TRUE & testval <= 0.05) {
                     grid.text("*", x, y, gp = gpar(fontsize = 14))
                   # NA values
                   } else if (is.na(mat_beta[i,j])) {
                     # grid.segments(
                     #   x0 = x + width/2, y0 = y + height/2,
                     #   x1 = x - width/2, y1 = y - height/2,
                     #   gp = gpar(col = "grey50", lwd = 1.1)
                     # )
                     grid.circle(x = x, y = y, r = unit(0.15, "cm"), gp = gpar(fill = "blue", col = NA))
                   }
                 })
  
  return(plt)
  
}


# -----------------------------------------------------------------------
# MAGMA ENRICHMENTS (PER CELLTYPE CONDITIONING)
# -----------------------------------------------------------------------

plot_heatmap_conditional_percelltype <- function(df, id = "BMI_2018", max_beta = NA, star = TRUE, title = NA,
                                                 named_colours = NULL) {
  
  # prep data
  dfx <- df |>
    filter(GWAS == id) |>
    mutate(CONDITIONED_ON = factor(CONDITIONED_ON, levels = names(named_colours$subcelltype))) |>
    arrange(CELLTYPE, CONDITIONED_ON)
  
  # set colours
  if (is.na(max_beta)) {
    max_beta <- ceiling(max(abs(dfx$BETA)))
  }
  col_fun <- circlize::colorRamp2(c(-max_beta, 0, max_beta), c("grey", "white", "red"))
  
  # create matrix (fdr and beta)
  mat_fdr <- dfx |>
    select(CELLTYPE, CONDITIONED_ON, FDR) |>
    tidyr::pivot_wider(names_from = CONDITIONED_ON, values_from = FDR, values_fill = Inf)
  mat_fdr <- mat_fdr[,-1]
  mat_fdr <- mat_fdr[, names(named_colours$subcelltype)]
  
  mat_beta <- dfx |>
    select(CELLTYPE, CONDITIONED_ON, BETA) |>
    tidyr::pivot_wider(names_from = CONDITIONED_ON, values_from = BETA, values_fill = list(value = NA_real_)) |>
    tibble::column_to_rownames(var = "CELLTYPE")  |>
    as.matrix()
  mat_beta <- mat_beta[, names(named_colours$subcelltype)]
  
  # create plot
  plt <- Heatmap(mat_beta, cluster_rows = F, cluster_columns = FALSE, col = col_fun,
                 width = dim(mat_beta)[2] * unit(0.7, "cm"),
                 height = dim(mat_beta)[1] * unit(0.8, "cm"),
                 row_names_side = "left",
                 na_col = "grey40",
                 column_title_gp = gpar(fontsize = 15),
                 column_names_gp = gpar(fontsize = 14),
                 row_names_gp = gpar(fontsize = 15),
                 column_title = ifelse(is.na(title), id, title),
                 # column_labels = c("-", new_conditioner$new_conditioner),
                 column_names_rot = 45,
                 heatmap_legend_param = list(
                   title = "BETA",
                   at = seq.int(-max_beta, max_beta, max_beta),
                   labels = seq.int(-max_beta, max_beta, max_beta),
                   labels_gp = gpar(fontsize = 12),
                   title_gp = gpar(fontsize = 12, fontface = "bold")
                 ), border = T,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   # signifant results (*, **)
                   testval <- mat_fdr[i, j]
                   if (star == TRUE & testval <= 0.01) {
                     grid.text("**", x, y, gp = gpar(fontsize = 14))
                   } else if (star == TRUE & testval <= 0.05) {
                     grid.text("*", x, y, gp = gpar(fontsize = 14))
                     # NA values
                   # } else if (is.na(mat_beta[i,j])) {
                   #   grid.segments(
                   #     x0 = x + width/2, y0 = y + height/2,
                   #     x1 = x - width/2, y1 = y - height/2,
                   #     gp = gpar(col = "grey50", lwd = 1.1)
                   #   )
                   }
                 })
  
  return(plt)
  
}



# -----------------------------------------------------------------------
# EWCE ENRICHMENT HEATMAPS
# -----------------------------------------------------------------------

# plot heatmap of EWCE enrichment results (sd from mean for colour, * for FDR)
plot_ewce_enrichment <- function(df, maxeffect = NA) {
  
  # set colours
  if (is.na(maxeffect)) {
    maxeffect <- ceiling(max(abs(df$sd_from_mean)))
  }
  col_fun <- circlize::colorRamp2(c(-maxeffect, 0, maxeffect), c("grey", "white", "red"))
  
  # matrix
  mat_effect <- df |>
    select(celltype, gene_set, sd_from_mean) |>
    tidyr::pivot_wider(names_from = celltype, values_from = sd_from_mean) |>
    tibble::column_to_rownames("gene_set")
  mat_effect <- mat_effect[, levels(df_ewce$celltype)] |>
    as.matrix()
  
  mat_fdr <- df |>
    select(celltype, gene_set, q) |>
    tidyr::pivot_wider(names_from = celltype, values_from = q) |>
    tibble::column_to_rownames("gene_set")
  mat_fdr <- mat_fdr[, levels(df_ewce$celltype)] |>
    as.matrix()
  
  p <- Heatmap(mat_effect, name = "ewce", cluster_rows = F, cluster_columns = F, col = col_fun,
               width = dim(mat_effect)[2] * unit(0.8, "cm"),
               height = dim(mat_effect)[1] * unit(0.8, "cm"),
               row_names_side = "left",
               column_names_rot = 45,
               column_title_gp = gpar(fontsize = 16),
               column_names_gp = gpar(fontsize = 12),
               column_title = df$gene_window[1],
               row_title = "T2D mechanistic cluster",
               row_title_gp = gpar(fontsize = 16),
               row_names_gp = gpar(fontsize = 12),
               heatmap_legend_param = list(title = "s.d", 
                                           title_gp = gpar(fontsize = 14, fontface = "bold"),
                                           at = c(-maxeffect, 0, maxeffect),
                                           labels_gp = gpar(fontsize = 14)),
               border = T,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 testval <- mat_fdr[i, j]
                 if (testval <= 0.01) {
                   grid.text("**", x, y, gp = gpar(fontsize = 14))
                 } else if (testval <= 0.05) {
                   grid.text("*", x, y, gp = gpar(fontsize = 14))
                 }
               })
  
  return(p)
}
