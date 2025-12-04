#### DESCRIPTION #######################################
# Purpose of script:
# Runs the individual metabolite MR analyses
#
# Author: Nick Sunderland
#
# Date Created: 2025-07-31
#
# Email: nicholas.sunderland@bristol.ac.uk

# #### SET INPUT #########################################
metab_outcome_mr_results <- snakemake@input[["metab_outcome_mr_results"]]
metab_outcome_pathway_mr <- snakemake@input[["metab_outcome_pathway_mr"]]
asp_hf_mr_plot     <- snakemake@output[["asp_hf_mr_plot"]]
asp_hf_mr_all_plot <- snakemake@output[["asp_hf_mr_all_plot"]]
# ########################################################

# requirements
library(data.table)
library(ggplot2)

metab_outcome_mr_results = "/Users/xx20081/git/bmi_metabolomics/output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz"
metab_outcome_pathway_mr = "/Users/xx20081/git/bmi_metabolomics/output/tables/mr_results/metabolite_pathway_outcome_mr_results.tsv.gz"
# read
outcome_dat <- fread(metab_outcome_mr_results)
pathway_dat <- fread(metab_outcome_pathway_mr)


# process
≠ <-

ggsave(asp_hf_mr_plot, p, width = 8, height = 5, dpi=300, bg="white")


# all MR methods
p <- ggplot(comb[!is.na(analysis)], aes(x = b, y = analysis, color = method, shape = p_sig)) +
  geom_vline(xintercept = 0, color="grey50", linetype="dashed") +
  geom_errorbarh(aes(xmin = b-1.96*b_se, xmax = b+1.96*b_se, group=method), height=0, position = position_dodge(width=.8)) +
  geom_point(position = position_dodge(width=.8), size = 2, fill = "white", stroke = 1.1) +
  scale_color_brewer(palette = "Set1", labels = ~ gsub("_", " ", .x)) +
  scale_shape_manual(values = c("P \u2265 0.05"=21, "P < 0.05"=19)) +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 9),
  ) +
  labs(x = "MR estimate - OR outcome per SD increase in exposure") +
  facet_wrap(~outcome, ncol=1, scales="free_y", labeller = labeller(outcome = function(x) {
    fcase(as.character(x)=="HFpEF", "Non-ischaemic HFpEF",
          as.character(x)=="HFrEF", "Non-ischaemic HFrEF",
          as.character(x)=="All-cause HF", "All-cause HF")
    }))
p

ggsave(asp_hf_mr_all_plot, p, width = 8, height = 8, dpi=300, bg="white")


