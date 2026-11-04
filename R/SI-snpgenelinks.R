# Supplementary Figures/Tables for Figure 4 (SNP-gene linking results)

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
  library(PRROC)
  library(argparse)
})

# parse arguments
parser <- ArgumentParser(description = "snMultiome figures")
parser$add_argument(
  "--config", required = TRUE, type = "character")

args <- parser$parse_args()

# config
source(args$config)

# -----------------------------------------------------------------------
# SET & CREATE DIRS/PATHS
# -----------------------------------------------------------------------

dir_in <- file.path(dir_data, "/15_manuscript_figures/data")
dir_out <- file.path(dir_scripts, "data/figures-supplementary/snp-gene-links")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# pgBoost links
df_links <- readRDS(paste0(dir_in, "/df_pgboost_links.rds"))
df_links <- df_links |>
  select(!grep("pgboost_", names(df_links), value = T))


# -----------------------------------------------------------------------
# AUPR
# -----------------------------------------------------------------------

# ref datasets to use
ref_ids <- grep("^in", names(df_links), value = T) |>
  setdiff("in_yengo_bmi")

# PR curves for each reference dataset
pr_list <- list()
pr_df_list <- list()

for (id in ref_ids) {
  pr <- pr.curve(scores.class0 = df_links$max_percentile_brain[df_links[[id]] == 1],
                 scores.class1 = df_links$max_percentile_brain[df_links[[id]] == 0],
                 curve = TRUE)
  baseline <- mean(df_links[[id]]) # proportion of links in this dataset
  pr_df <- as.data.frame(pr$curve) |>
    rename(recall = V1, precision = V2, threshold = V3) |>
    mutate(enrichment = precision / baseline)
  pr_list[[id]] <- pr
  pr_df_list[[id]] <- pr_df
}

# Numbers of links per reference
table_links <- df_links |>
  mutate(in_candidate_links = 1) |>
  tidyr::pivot_longer(cols = starts_with("in_"), names_to = "dataset", values_to = "in_dataset") |>
  mutate(dataset = gsub("in_", "", dataset)) |>
  group_by(dataset) |>
  summarise(n_links = sum(in_dataset), n_total = n(), pc_links = 100 * n_links / n_total) |>
  select(-n_total) |>
  arrange(desc(pc_links))

write.csv(table_links, file.path(dir_out, "table_validation_sets_comparison.csv"), row.names = F)

# -----------------------------------------------------------------------
# FIGURE SX: AUPR CURVES
# -----------------------------------------------------------------------

# Study names
ref_names <- read.csv(file.path(dir_scripts, "data/metadata/snp-gene-link-references.csv"))

plot_names <- setNames(ref_names$name, ref_names$id)
names(pr_list) <- gsub("in_", "", names(pr_list))
names(pr_df_list) <- gsub("in_", "", names(pr_df_list))

# Panel per ref id
plot_list <- sapply(names(plot_names), function(id) {
  p <- ggplot(pr_df_list[[id]], aes(x = recall, y = enrichment)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    labs(x = "Recall", y = "Enrichment", title = plot_names[id]) +
    theme_classic() +
    theme(plot.title = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)) +
    scale_x_continuous(breaks = c(0, 0.5, 1))
})

# Figure
pdf(file.path(dir_out, "fig_aupr_curves.pdf"), width = 6.5, height = 4.5)
wrap_plots(plot_list, ncol = 3, axis_titles = "collect")
invisible(dev.off())


# -----------------------------------------------------------------------
# TABLE SX: AUPR DATA
# -----------------------------------------------------------------------

df_aupr <- ref_names |>
  left_join(table_links, join_by(id == dataset)) |>
  mutate(baseline = pc_links / 100) |>
  select(!pc_links)

df_aupr$aupr <- sapply(df_aupr$id, function(i) {pr_list[[i]]$auc.integral})

write.csv(df_aupr, file.path(dir_out, "table_aupr_results.csv"), row.names = F)
