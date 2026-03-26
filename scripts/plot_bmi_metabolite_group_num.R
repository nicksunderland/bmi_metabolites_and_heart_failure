#### DESCRIPTION #######################################
# Purpose of script:
# Plot the trial consistent metabolites by pathway group
#
# Author: Nick Sunderland
#
# Date Created: 2026-03-17
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
replicating    <- snakemake@input[["replicating"]]
map            <- snakemake@input[["map"]]
rep_by_grp     <- snakemake@input[["rep_by_grp"]]
fig_rep_by_grp <- snakemake@output[["fig_rep_by_grp"]]
########################################################

# testing
if (FALSE) {
  replicating    <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/replication/metabolite_replication_bbs_direct_mr.tsv")
  map            <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "scripts/gwas_metab_name_map.xlsx")
  rep_by_grp     <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/replication/metabolite_replication_bbs_direct_mr_by_group.tsv")
  fig_rep_by_grp <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/figures/replication/metabolite_replication_bbs_direct_mr_by_group.png")
}

# required packages
library(data.table)
library(readxl)
library(ggplot2)

# read in data
res <- fread(replicating)
map <- read_xlsx(map) |> as.data.table()

# clean
res <- res[replicates_trials_mr == "Trial & MR consistent"]
res[map, pathway := i.pathway_group, on = "label"]
res[, pathway := factor(pathway, levels = res[, .N, by = pathway][order(-N), pathway])]

# write out the numbers by group
fwrite(res, rep_by_grp, sep="\t")

# plot the numbers by group
p <- ggplot(res, aes(y = pathway, fill = pathway)) +
  geom_bar() +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    hjust = -0.3,
    size = 3.5
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.75, 0.75)
  ) +
  labs(
    fill = "Metabolite group",
    x = "Count",
    y = NULL
  )
p

# save the figure
ggsave(
  fig_rep_by_grp,
  p,
  dpi=300,
  width = 7,
  height = 7,
  bg = "white"
)

#
# library(ggrepel)
#
# # res[, bbs_bmi_decr_scaled := -as.numeric(bmi_estimate_bbs)]
# # setcolorder(res, "bbs_bmi_decr_scaled")
#
# p <- ggplot(res, aes(x = -as.numeric(bmi_estimate_bbs), y = as.numeric(bmi_estimate_mr))) +
#   geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
#   geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
#   geom_point() +
#   geom_point(data = res[sign(as.numeric(bmi_estimate_bbs)) != sign(as.numeric(bmi_estimate_mr)) &
#                           bmi_p.value_mr < 0.05, ], color="red") +
#   geom_label_repel(data = res[sign(as.numeric(bmi_estimate_bbs)) != sign(as.numeric(bmi_estimate_mr)) &
#                                 bmi_p.value_mr < 0.05, ],
#                    aes(label = label),
#                    size = 3,
#                    max.overlaps = Inf) +
#   geom_smooth(method="lm", se=F) +
#   theme_classic() +
#   labs(x = "Weight loss intervention (\u0394metabolite per SD BMI decrease)",
#        y = "Lifetime BMI exposure (\u0394metabolite per SD BMI increase)")
#
# ggsave(
#   "/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/.results/figures/replication/simple_corr_bmi_surgery.png",
#   p,
#   dpi=300,
#   width = 6,
#   height = 5.5,
#   bg = "white"
# )
#


