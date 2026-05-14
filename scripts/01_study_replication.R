#### DESCRIPTION #######################################
# Compare weight-loss trial association statistics between
# BBS and DiRECT studies.
# Produces: per-study volcano plots, within-study platform
# correlation (MS vs NMR), and BBS vs DiRECT replication
# scatter / Bland-Altman / summary table.
########################################################

#### SET INPUT #########################################
assoc           <- snakemake@input[["assoc"]]
name_map_file   <- snakemake@input[["name_map"]]
assoc_table     <- snakemake@output[["assoc_table"]]
vol_bbs         <- snakemake@output[["vol_bbs"]]
vol_direct      <- snakemake@output[["vol_direct"]]
plat_corr       <- snakemake@output[["plat_corr"]]
plat_bland_alt  <- snakemake@output[["plat_bland_alt"]]
corr_plot       <- snakemake@output[["corr_plot"]]
study_bland_alt <- snakemake@output[["study_bland_alt"]]
rep_plot        <- snakemake@output[["rep_plot"]]
rep_table       <- snakemake@output[["rep_table"]]
consistent_ids       <- snakemake@output[["consistent_ids"]]
direct_est_scatter   <- snakemake@output[["direct_est_scatter"]]
direct_est_bland_alt <- snakemake@output[["direct_est_bland_alt"]]
direct_se_scatter    <- snakemake@output[["direct_se_scatter"]]
log_file             <- snakemake@log[["log"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
repo_dir <- Sys.getenv("HF_METABOLITE_REPO2")
assoc <- c(
  file.path(repo_dir, "output/tables/linear_mixed_associations/linear_mixed_associations_raw_z_bbs_nightingale.tsv"),
  file.path(repo_dir, "output/tables/linear_mixed_associations/linear_mixed_associations_raw_z_direct_nightingale.tsv"),
  file.path(repo_dir, "output/tables/linear_mixed_associations/linear_mixed_associations_rnt_bbs_metabolon.tsv"),
  file.path(repo_dir, "output/tables/linear_mixed_associations/linear_mixed_associations_rnt_direct_metabolon.tsv")
)
assoc_table     <- file.path(repo_dir, "output/tables/linear_mixed_associations/combined_trial_associations.tsv")
vol_bbs         <- file.path(repo_dir, "output/figures/associations/bbs_linear_mixed_volcanos.png")
vol_direct      <- file.path(repo_dir, "output/figures/associations/direct_linear_mixed_volcanos.png")
plat_corr       <- file.path(repo_dir, "output/figures/replication/platform_measurement_correlation.png")
plat_bland_alt  <- file.path(repo_dir, "output/figures/replication/platform_measurement_bland_alterman.png")
corr_plot       <- file.path(repo_dir, "output/figures/replication/study_comparison_correlation.png")
study_bland_alt <- file.path(repo_dir, "output/figures/replication/study_comparison_bland_alterman.png")
rep_plot        <- file.path(repo_dir, "output/figures/replication/study_comparison_replication.png")
rep_table       <- file.path(repo_dir, "output/tables/replication/study_replication.tsv")
name_map_file   <- file.path(repo_dir, "scripts/gwas_metab_name_map.xlsx")
consistent_ids       <- file.path(repo_dir, "output/tables/replication/consistent_gwas_ids.txt")
direct_est_scatter   <- file.path(repo_dir, "output/figures/replication/direct_estimate_comparison_scatter.png")
direct_est_bland_alt <- file.path(repo_dir, "output/figures/replication/direct_estimate_comparison_bland_altman.png")
direct_se_scatter    <- file.path(repo_dir, "output/figures/replication/direct_se_comparison_scatter.png")
log_file             <- file.path(repo_dir, "output/logs/study_replication.log")
# ─────────────────────────────────────────────────────────────────────────────
}

library(ggplot2)
library(ggrepel)
library(patchwork)
library(mcr)
library(readxl)
suppressPackageStartupMessages(library(data.table))

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
on.exit(close(log_con))
tee <- function(...) { cat(...); cat(..., file = log_con, append = TRUE) }


#############
# WEIGHT LOSS TRIAL ASSOCIATIONS
#############

res <- lapply(assoc, fread) |> rbindlist(use.names = TRUE, fill = TRUE)

res[grepl("^olink_", feature_id), platform := "olink"]
total_metabs_assessed <- res[, .(n = uniqueN(label)), by = .(main_study, platform)]
rbind(
  total_metabs_assessed[, .(platform = "all", n = sum(n)), by = "main_study"],
  total_metabs_assessed
)

res <- res[term == "timepointend" & model_name == "linear_mixed_time"]

# p value adjustments made within metabolite platform types x study
# - but should be within study (not also within measurement platform) - recompute
res[, p.value_fdr := p.adjust(p.value, method = "fdr", n = .N), by = .(main_study)]
res[, p.value_holm := p.adjust(p.value, method = "holm", n = .N), by = .(main_study)]
res[, p.value_interaction_fdr := p.adjust(p.value_interaction, method = "fdr", n = .N), by = .(main_study)]
res[, p.value_interaction_holm := p.adjust(p.value_interaction, method = "holm", n = .N), by = .(main_study)]


res[, `:=`(main_study    = as.factor(main_study),
           platform      = as.factor(platform),
           pathway_group = factor(pathway_group, levels = c(
             "Lipid metabolism", "Lipoproteins", "Amino acid metabolism",
             "Carbohydrate/energy metabolism", "Nucleotide metabolism",
             "Cofactors and vitamins", "Inflammation", "Fluid balance",
             "Xenobiotics/other", "Unknown")))]

n_top <- 10
res   <- res[order(p.value)]
res[, top_label := ifelse(seq_len(.N) < n_top, label, NA_character_), by = .(main_study, platform)]


# wide table: both trials side by side
assoc_wide <- merge(
  res[main_study == "bbs", .(
    label, platform, pathway_group,
    bbs_estimate = estimate,
    bbs_se       = std.error,
    bbs_p        = p.value,
    bbs_p_fdr    = p.value_fdr,
    bbs_n        = n,
    bbs_model    = formula
  )],
  res[main_study == "direct", .(
    label, platform,
    # intervention only model estimates ===
    direct_estimate           = estimate,
    direct_se                 = std.error,
    direct_p                  = p.value,
    direct_p_fdr              = p.value_fdr,
    direct_n                  = n,
    direct_model              = formula,
    # full interaction model ===
    # timepoint main effect (change in metab in the controls)
    direct_estimate_time        = estimate_time,
    direct_se_time              = se_time,
    direct_p_time               = p.value_time,
    direct_p_time_fdr           = p.value_time_fdr,
    # time-allocation interaction effect (additional change in the intervention group)
    direct_estimate_interaction = estimate_interaction,
    direct_se_interaction       = se_interaction,
    direct_p_interaction        = p.value_interaction,
    direct_p_interaction_fdr    = p.value_interaction_fdr,
    # full model info
    direct_full_n               = n_full,
    direct_full_model           = formula_full
  )],
  by = c("label", "platform"), all = TRUE
)
setorder(assoc_wide, pathway_group, label, platform)

dir.create(dirname(assoc_table), recursive = TRUE, showWarnings = FALSE)
fwrite(assoc_wide, assoc_table, sep = "\t")


# print summary counts
tee(sprintf(
  "DiRECT significant associations:
  total=%i; nominal p<0.05: n=%i (%.1f%%); FDR p<0.05: n=%i (%.1f%%); Holm p<0.05: n=%i (%.1f%%)
  Nominal: increasing=%i; decreasing=%i
  FDR:     increasing=%i; decreasing=%i\n",
  total_metabs_assessed[main_study == "direct", sum(n)],
  res[main_study == "direct" & !is.na(estimate), sum(p.value_interaction < 0.05)],
  res[main_study == "direct" & !is.na(estimate), 100 * sum(p.value_interaction < 0.05) / total_metabs_assessed[main_study == "direct", sum(n)]],
  res[main_study == "direct" & !is.na(estimate), sum(p.value_interaction_fdr < 0.05)],
  res[main_study == "direct" & !is.na(estimate), 100 * sum(p.value_interaction_fdr < 0.05) / total_metabs_assessed[main_study == "direct", sum(n)]],
  res[main_study == "direct" & !is.na(estimate), sum(p.value_interaction_holm < 0.05)],
  res[main_study == "direct" & !is.na(estimate), 100 * sum(p.value_interaction_holm < 0.05) / total_metabs_assessed[main_study == "direct", sum(n)]],
  res[main_study == "direct" & !is.na(estimate) & p.value_interaction < 0.05, sum(sign(estimate) == 1)],
  res[main_study == "direct" & !is.na(estimate) & p.value_interaction < 0.05, sum(sign(estimate) == -1)],
  res[main_study == "direct" & !is.na(estimate) & p.value_interaction_fdr < 0.05, sum(sign(estimate) == 1)],
  res[main_study == "direct" & !is.na(estimate) & p.value_interaction_fdr < 0.05, sum(sign(estimate) == -1)]
))

tee(sprintf(
  "BBS significant associations:
  total=%i; nominal p<0.05: n=%i (%.1f%%); FDR p<0.05: n=%i (%.1f%%); Holm p<0.05: n=%i (%.1f%%)
  Nominal: increasing=%i; decreasing=%i
  FDR:     increasing=%i; decreasing=%i\n",
  total_metabs_assessed[main_study == "bbs", sum(n)],
  res[main_study == "bbs" & !is.na(estimate), sum(p.value < 0.05)],
  res[main_study == "bbs" & !is.na(estimate), 100 * sum(p.value < 0.05) / total_metabs_assessed[main_study == "bbs", sum(n)]],
  res[main_study == "bbs" & !is.na(estimate), sum(p.value_fdr < 0.05)],
  res[main_study == "bbs" & !is.na(estimate), 100 * sum(p.value_fdr < 0.05) / total_metabs_assessed[main_study == "bbs", sum(n)]],
  res[main_study == "bbs" & !is.na(estimate), sum(p.value_holm < 0.05)],
  res[main_study == "bbs" & !is.na(estimate), 100 * sum(p.value_holm < 0.05) / total_metabs_assessed[main_study == "bbs", sum(n)]],
  res[main_study == "bbs" & !is.na(estimate) & p.value < 0.05, sum(sign(estimate) == 1)],
  res[main_study == "bbs" & !is.na(estimate) & p.value < 0.05, sum(sign(estimate) == -1)],
  res[main_study == "bbs" & !is.na(estimate) & p.value_fdr < 0.05, sum(sign(estimate) == 1)],
  res[main_study == "bbs" & !is.na(estimate) & p.value_fdr < 0.05, sum(sign(estimate) == -1)]
))


# volcano plots
volcanos <- list(bbs = vol_bbs, direct = vol_direct)
for (s in names(volcanos)) {
  df      <- res[main_study == s]
  df$yval <- if (s == "bbs") -log10(df$p.value) else -log10(df$p.value) #_interaction

  p <- ggplot(df, aes(x = estimate, y = yval)) +
    geom_point(aes(fill = pathway_group), shape = 21, color = "black",
               size = 2, stroke = 0.3, alpha = 0.7) +
    geom_text_repel(aes(label = top_label), size = 2.5, color = "black") +
    theme_light(base_size = 12) +
    labs(
      y    = if (s == "bbs") expression(-log[10](p.value)) else expression(-log[10](p.value~interaction)),
      x    = expression(beta[timepoint]),
      fill = "Pathway"
    ) +
    theme(panel.grid = element_blank(), plot.margin = margin(t = 30, r = 5, b = 5, l = 5)) +
    facet_wrap(~ platform, ncol = 1, scales = "free",
               labeller = labeller(platform = function(x) ifelse(x == "nightingale", "NMR", "MS")))

  dir.create(dirname(volcanos[[s]]), recursive = TRUE, showWarnings = FALSE)
  ggsave(volcanos[[s]], p, width = 8, height = 8, dpi = 300, bg = "white")
}


#############
# PLATFORM CORRELATION (MS vs NMR within study)
#############
in_both_platforms <- res[, count_name := .N, by = c("label", "main_study")][count_name == 2]
in_both_platforms <- dcast(in_both_platforms, label + main_study ~ platform,
                            value.var = c("estimate", "std.error"))

deming_dat <- in_both_platforms[, {
  error_ratio <- mean(std.error_metabolon^2,   na.rm = TRUE) /
                 mean(std.error_nightingale^2, na.rm = TRUE)
  fit        <- mcr::mcreg(estimate_metabolon, estimate_nightingale,
                            method.reg = "Deming", error.ratio = error_ratio)
  cor_ci     <- confintr::ci_cor(estimate_metabolon, estimate_nightingale)
  pearson    <- cor_ci[["estimate"]]
  pearson_lb <- cor_ci[["interval"]][[1]]
  pearson_ub <- cor_ci[["interval"]][[2]]
  .(slope     = coef(fit)[2],
    intercept = coef(fit)[1],
    lab_x     = -1,
    lab_y     = ifelse(.BY[["main_study"]] == "bbs", 0.5, 0.4),
    label     = sprintf("Pearson's R: %.2f (%.2f-%.2f) (%s)",
                        pearson, pearson_lb, pearson_ub, toupper(.BY[["main_study"]])))
}, by = "main_study"]

p <- ggplot(in_both_platforms, aes(x = estimate_metabolon, y = estimate_nightingale, color = main_study)) +
  geom_errorbar(aes(ymin = estimate_nightingale - 1.96 * std.error_nightingale,
                    ymax = estimate_nightingale + 1.96 * std.error_nightingale),
                color = "lightgray", width = 0) +
  geom_errorbar(aes(xmin = estimate_metabolon - 1.96 * std.error_metabolon,
                    xmax = estimate_metabolon + 1.96 * std.error_metabolon),
                orientation = "y",
                color = "lightgray", height = 0) +
  geom_abline(data = deming_dat, aes(intercept = intercept, slope = slope, color = main_study)) +
  geom_text_repel(aes(label = label), color = "black") +
  geom_point() +
  geom_text(data = deming_dat, aes(x = lab_x, y = lab_y, label = label),
            fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c(bbs = "#8DA0CB", direct = "#FC8D62"),
                     labels = function(x) toupper(x)) +
  theme_classic() +
  labs(x = "MS estimate", y = "NMR estimate", color = "Study")
p

dir.create(dirname(plat_corr), recursive = TRUE, showWarnings = FALSE)
ggsave(plat_corr, p, width = 8, height = 6, dpi = 300, bg = "white")


# Bland-Altman (platform)
in_both_platforms[, `:=`(
  mean_estimate = (estimate_metabolon + estimate_nightingale) / 2,
  diff_estimate =  estimate_metabolon - estimate_nightingale
)]
ba_plat <- in_both_platforms[, .(`  `                       = 0,
                                  Bias                       = mean(diff_estimate),
                                  `Upper limit of agreement` = mean(diff_estimate) + 1.96 * sd(diff_estimate),
                                  `Lower limit of agreement` = mean(diff_estimate) - 1.96 * sd(diff_estimate))]
ba_plat <- melt(ba_plat)
ba_plat[, linetype := ifelse(variable == "  ", "solid", "dashed")]

p <- ggplot(in_both_platforms, aes(x = mean_estimate, y = diff_estimate, fill = main_study)) +
  geom_hline(data = ba_plat, aes(yintercept = value, linetype = linetype)) +
  geom_point(shape = 21, col = "black", size = 3) +
  geom_text(data = ba_plat,
            aes(label = variable,
                x     = max(in_both_platforms$mean_estimate, na.rm = TRUE) + 0.02,
                y     = value + 0.02),
            hjust = 1, color = "gray", inherit.aes = FALSE, size = 3) +
  geom_text_repel(aes(label = label), color = "black", size = 3) +
  scale_fill_manual(values = c(bbs = "#8DA0CB", direct = "#FC8D62"),
                    labels = function(x) toupper(x)) +
  scale_linetype_identity() +
  theme_light() +
  theme(legend.position = "top") +
  labs(y = "Difference", x = "Mean", fill = "Study")
p

ggsave(plat_bland_alt, p, width = 8, height = 6, dpi = 300, bg = "white")


#############
# DiRECT: INT-ONLY vs COMBINED ESTIMATE VALIDATION
# Demonstrates that the intervention-arm-only model (estimate) and the
# combined full-model estimate (time + interaction) are consistent
#############
direct_res <- res[main_study == "direct" & !is.na(estimate) & !is.na(estimate_combined)]

deming_dc <- direct_res[!is.na(std.error) & !is.na(se_combined), {
  er  <- mean(std.error^2) / mean(se_combined^2)
  fit <- mcr::mcreg(estimate, estimate_combined, method.reg = "Deming", error.ratio = er)
  cc  <- confintr::ci_cor(estimate, estimate_combined)
  .(slope     = coef(fit)[2],
    intercept = coef(fit)[1],
    lab       = sprintf("r = %.2f (%.2f–%.2f)", cc$estimate, cc$interval[1], cc$interval[2]))
}, by = "platform"]

direct_res[deming_dc, `:=`(dc_slope = i.slope, dc_int = i.intercept), on = "platform"]
direct_res[, dc_orth := abs(dc_slope * estimate - estimate_combined + dc_int) / sqrt(dc_slope^2 + 1)]
direct_res[, dc_flag := fifelse(frank(-dc_orth, ties.method = "first") <= 8,
                                label, NA_character_), by = "platform"]

# 1. Scatter: int-only vs combined
p_dc_scatter <- ggplot(direct_res, aes(x = estimate, y = estimate_combined, colour = platform)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linewidth = 0.35, linetype = "dashed") +
  geom_abline(data = deming_dc, aes(intercept = intercept, slope = slope, colour = platform),
              linewidth = 0.7) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text(data = deming_dc, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.5, fontface = "bold", size = 3, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  guides(colour = "none") +
  facet_wrap(~ platform, ncol = 1,
             labeller = labeller(platform = function(x) ifelse(x == "nightingale", "NMR", "MS"))) +
  labs(x = expression(beta ~ "(DiRECT, intervention-arm only model)"),
       y = expression(beta ~ "(DiRECT, time + interaction from full model)")) +
  theme_light(base_size = 12)

dir.create(dirname(direct_est_scatter), recursive = TRUE, showWarnings = FALSE)
ggsave(direct_est_scatter, p_dc_scatter, width = 7, height = 10, dpi = 300, bg = "white")


# 2. Bland-Altman: int-only vs combined
direct_res[, `:=`(dc_ba_mean = (estimate + estimate_combined) / 2,
                   dc_ba_diff =  estimate - estimate_combined)]

dc_ba_lines <- direct_res[, .(`  `                       = 0,
                                Bias                       = mean(dc_ba_diff, na.rm = TRUE),
                                `Upper limit of agreement` = mean(dc_ba_diff, na.rm = TRUE) + 1.96 * sd(dc_ba_diff, na.rm = TRUE),
                                `Lower limit of agreement` = mean(dc_ba_diff, na.rm = TRUE) - 1.96 * sd(dc_ba_diff, na.rm = TRUE))]
dc_ba_lines <- melt(dc_ba_lines)
dc_ba_lines[, linetype := ifelse(variable == "  ", "solid", "dashed")]

p_dc_ba <- ggplot(direct_res, aes(x = dc_ba_mean, y = dc_ba_diff, fill = platform)) +
  geom_hline(data = dc_ba_lines, aes(yintercept = value, linetype = linetype)) +
  geom_point(shape = 21, stroke = 0, size = 2, alpha = 0.7) +
  geom_text(data = dc_ba_lines,
            aes(label = variable,
                x     = max(direct_res$dc_ba_mean, na.rm = TRUE),
                y     = value),
            hjust = 1, vjust = -0.5, colour = "gray", inherit.aes = FALSE, size = 3) +
  scale_fill_brewer(palette = "Set1",
                    labels = function(x) ifelse(x == "nightingale", "NMR", "MS")) +
  scale_linetype_identity() +
  theme_light() +
  theme(legend.position = "top") +
  labs(y = "Difference (int-only − combined)", x = "Mean", fill = "Platform")

ggsave(direct_est_bland_alt, p_dc_ba, width = 10, height = 7, dpi = 300, bg = "white")


# 3. SE comparison: int-only vs combined (approx)
p_dc_se <- ggplot(direct_res[!is.na(std.error) & !is.na(se_combined)],
                  aes(x = std.error, y = se_combined, colour = platform)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linewidth = 0.35, linetype = "dashed") +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_brewer(palette = "Set1",
                      labels = function(x) ifelse(x == "nightingale", "NMR", "MS")) +
  guides(colour = "none") +
  facet_wrap(~ platform, ncol = 1,
             labeller = labeller(platform = function(x) ifelse(x == "nightingale", "NMR", "MS"))) +
  labs(x = "SE (intervention-arm only model)",
       y = "SE (combined, approx.)") +
  theme_light(base_size = 12)

ggsave(direct_se_scatter, p_dc_se, width = 6, height = 9, dpi = 300, bg = "white")


#############
# BBS vs DiRECT REPLICATION
#############
replication <- dcast(res, label + platform + pathway_group ~ main_study,
                     value.var = c("estimate", "std.error", "p.value", "p.value_fdr",
                                   "p.value_interaction", "p.value_interaction_fdr",
                                   "p.value_holm", "n"))

replication[, `:=`(
  same_direction     = sign(estimate_bbs) == sign(estimate_direct),
  both_significant   = p.value_fdr_bbs < 0.05 & p.value_interaction_fdr_direct < 0.05,
  interaction_backed = !is.na(p.value_interaction_fdr_direct) & p.value_interaction_fdr_direct < 0.05
)]

replication[, replicates := factor(
  fcase(same_direction & both_significant & interaction_backed, "Trial consistent (interaction)",
        same_direction & both_significant,                       "Trial consistent",
        both_significant,                                        "Trial FDR P < 0.05 only",
        same_direction,                                          "Trial direction only",
        default = "None"),
  levels = c("Trial consistent (interaction)", "Trial consistent",
             "Trial FDR P < 0.05 only", "Trial direction only", "None")
)]

setkeyv(replication, c("replicates", "label"))

tee(sprintf(
  "Total features=%i (BBS=%i; DiRECT=%i), Overlap n=%i (%.1f%%), concordant=%i (%.1f%%), discordant=%i (%.1f%%)\n",
  replication[, uniqueN(interaction(label, platform))],
  replication[!is.na(estimate_bbs),    uniqueN(interaction(label, platform))],
  replication[!is.na(estimate_direct), uniqueN(interaction(label, platform))],
  replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
  100 * replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)] / replication[, uniqueN(label)],
  replication[, sum(same_direction,  na.rm = TRUE)],
  100 * replication[, sum(same_direction,  na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
  replication[, sum(!same_direction, na.rm = TRUE)],
  100 * replication[, sum(!same_direction, na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)]
))

consistent_levels <- c("Trial consistent (interaction)", "Trial consistent")
tee(sprintf(
  "Replication: sign=%i (%.1f%%), sig=%i (%.1f%%), both=%i (%.1f%%)
  Trial consistent +ve=%i (%.1f%%); -ve=%i (%.1f%%)
  Of Trial consistent: interaction-backed n=%i (%.1f%%)\n",
  replication[, sum(same_direction,   na.rm = TRUE)],
  100 * replication[, sum(same_direction,   na.rm = TRUE)] / replication[, sum(!is.na(same_direction))],
  replication[, sum(both_significant, na.rm = TRUE)],
  100 * replication[, sum(both_significant, na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct))],
  replication[replicates %in% consistent_levels, .N],
  100 * replication[replicates %in% consistent_levels, .N] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
  replication[replicates %in% consistent_levels & sign(estimate_bbs) ==  1, .N],
  100 * replication[replicates %in% consistent_levels & sign(estimate_bbs) ==  1, .N] / replication[replicates %in% consistent_levels, .N],
  replication[replicates %in% consistent_levels & sign(estimate_bbs) == -1, .N],
  100 * replication[replicates %in% consistent_levels & sign(estimate_bbs) == -1, .N] / replication[replicates %in% consistent_levels, .N],
  replication[replicates == "Trial consistent (interaction)", .N],
  100 * replication[replicates == "Trial consistent (interaction)", .N] / replication[replicates %in% consistent_levels, .N]
))

replication[, uniqueN(interaction(label, platform)), by = "replicates"]

# Pathway group breakdown of replicating metabolites
n_consistent <- replication[replicates %in% consistent_levels, .N]
pg_counts <- replication[replicates %in% consistent_levels,
                          .N, by = pathway_group][order(-N)]
tee("Replicating metabolites by pathway group (% of replicating):\n")
for (i in seq_len(nrow(pg_counts))) {
  tee(sprintf("  %-35s n=%i (%.1f%%)\n",
              as.character(pg_counts$pathway_group[i]),
              pg_counts$N[i],
              100 * pg_counts$N[i] / n_consistent))
}


# BBS vs DiRECT correlation
corr_dat <- dcast(res[model_name == "linear_mixed_time"][abs(estimate) > 2, estimate := NA_real_],
                  label + platform ~ main_study, value.var = c("estimate", "std.error"))
corr_dat <- corr_dat[complete.cases(corr_dat)]

deming_corr <- corr_dat[, {
  error_ratio <- mean(std.error_bbs^2,    na.rm = TRUE) /
                 mean(std.error_direct^2, na.rm = TRUE)
  fit        <- mcr::mcreg(estimate_bbs, estimate_direct,
                            method.reg = "Deming", error.ratio = error_ratio)
  cor_ci     <- confintr::ci_cor(estimate_bbs, estimate_direct)
  pearson    <- cor_ci[["estimate"]]
  pearson_lb <- cor_ci[["interval"]][[1]]
  pearson_ub <- cor_ci[["interval"]][[2]]
  .(slope     = coef(fit)[2],
    intercept = coef(fit)[1],
    lab_x     = -1.2,
    lab_y     = 1.15,
    label     = sprintf("Pearson's R: %.2f (%.2f-%.2f)", pearson, pearson_lb, pearson_ub))
}, by = "platform"]

corr_dat[deming_corr, `:=`(slope = i.slope, intercept = i.intercept), on = "platform"]
corr_dat[, orth_resid     := abs(slope * estimate_bbs - estimate_direct + intercept) / sqrt(slope^2 + 1)]
corr_dat[, labels_to_plot := fifelse(frank(-orth_resid, ties.method = "first") <= 10,
                                     label, NA_character_), by = "platform"]

p <- ggplot(corr_dat, aes(x = estimate_bbs, y = estimate_direct, color = platform)) +
  geom_point() +
  geom_abline(data = deming_corr, aes(intercept = intercept, slope = slope, color = platform)) +
  geom_text_repel(aes(label = labels_to_plot), color = "black", size = 2.5, na.rm = TRUE) +
  geom_text(data = deming_corr, aes(x = lab_x, y = lab_y, label = label),
            fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  theme_light(base_size = 14) +
  labs(x = expression(beta[timepoint~BBS]), y = expression(beta[timepoint~DiRECT])) +
  guides(color = "none") +
  facet_wrap(~ platform, ncol = 1,
             labeller = labeller(platform = function(x) ifelse(x == "nightingale", "NMR", "MS")))
p

dir.create(dirname(corr_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(corr_plot, p, width = 7, height = 10, dpi = 300)


# Bland-Altman (study)
corr_dat[, `:=`(
  mean_estimate = (estimate_bbs + estimate_direct) / 2,
  diff_estimate =  estimate_bbs - estimate_direct
)]
ba_study <- corr_dat[, .(`  `                       = 0,
                          Bias                       = mean(diff_estimate),
                          `Upper limit of agreement` = mean(diff_estimate) + 1.96 * sd(diff_estimate),
                          `Lower limit of agreement` = mean(diff_estimate) - 1.96 * sd(diff_estimate))]
ba_study <- melt(ba_study)
ba_study[, linetype := ifelse(variable == "  ", "solid", "dashed")]

p <- ggplot(corr_dat, aes(x = mean_estimate, y = diff_estimate, fill = platform)) +
  geom_hline(data = ba_study, aes(yintercept = value, linetype = linetype)) +
  geom_point(shape = 21, stroke = 0, size = 2, alpha = 0.7) +
  geom_text(data = ba_study,
            aes(label = variable,
                x     = max(corr_dat$mean_estimate, na.rm = TRUE),
                y     = value),
            hjust = 1, vjust = -0.5, color = "gray", inherit.aes = FALSE, size = 3) +
  geom_text_repel(
    data = corr_dat[diff_estimate > ba_study[variable == "Upper limit of agreement", value] |
                    diff_estimate < ba_study[variable == "Lower limit of agreement", value]],
    aes(label = label), color = "black", size = 3) +
  scale_fill_brewer(palette = "Set1",
                    labels = function(x) ifelse(x == "nightingale", "NMR", "MS")) +
  scale_linetype_identity() +
  theme_light() +
  theme(legend.position = "top") +
  labs(y = "Difference", x = "Mean", fill = "Platform")
p

ggsave(study_bland_alt, p, width = 10, height = 7, dpi = 300, bg = "white")


# replication scatter
replication[, pathway_broad := fcase(
  pathway_group == "Lipid metabolism", "Lipid metabolism",
  pathway_group == "Lipoproteins",     "Lipoproteins",
  pathway_group == "Amino acid metabolism", "Amino acids",
  default = "Other"
)]
replication[, pathway_broad := factor(pathway_broad,
  levels = c("Lipid metabolism", "Lipoproteins", "Amino acids", "Other"))]

p <- ggplot(replication, aes(x = estimate_bbs, y = estimate_direct)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(
    data   = replication[replicates %in% c("None", "Trial direction only")],
    colour = "grey78", size = 0.8, alpha = 0.45
  ) +
  geom_point(
    data  = replication[replicates %in% c("Trial consistent (interaction)", "Trial consistent", "Trial FDR P < 0.05 only")],
    aes(fill = replicates), shape = 21, colour = "grey25", stroke = 0.3, size = 2.5
  ) +
  geom_text_repel(
    data = {
      consistent <- replication[replicates %in% consistent_levels]
      consistent[, effect_score := sqrt(abs(estimate_bbs) * abs(estimate_direct))]
      opposing   <- replication[replicates == "Trial FDR P < 0.05 only"]
      rbind(
        consistent[consistent[, frank(-effect_score, ties.method = "first") <= 10, by = pathway_broad]$V1][, effect_score := NULL],
        opposing
      )
    },
    aes(label          = label),
    size               = 2.5,
    colour             = "grey15",
    max.overlaps       = 20,
    box.padding        = 0.4,
    min.segment.length = 0.2,
    force = 20,
    segment.colour     = "grey55",
    segment.size       = 0.3,
    bg.color = alpha("white", 0.6),
    bg.r     = 0.08
  ) +
  scale_fill_manual(
    values = c("Trial consistent (interaction)" = "#B2182B",
               "Trial consistent"               = "#D6604D",
               "Trial FDR P < 0.05 only"       = "#2166AC"),
    labels = c("Trial consistent (interaction)" = "FDR P < 0.05, consistent + DiRECT interaction",
               "Trial consistent"               = "FDR P < 0.05, consistent (int-only)",
               "Trial FDR P < 0.05 only"       = "FDR P < 0.05, opposing direction")
  ) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  facet_wrap(~ pathway_broad, ncol = 2) +
  labs(
    x    = expression(beta[timepoint]~"(BBS)"),
    y    = expression(beta[timepoint]~"(DiRECT)"),
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    legend.key.size  = unit(0.8, "lines"),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(8, 8, 8, 8)
  )
p

dir.create(dirname(rep_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(rep_plot, p, width = 10, height = 10, dpi = 300, bg = "white")


# save replication table
dir.create(dirname(rep_table), recursive = TRUE, showWarnings = FALSE)
fwrite(replication[, .(
  label,
  platform,
  replicates,
  interaction_backed,
  estimate_bbs,
  estimate_direct,
  std.error_bbs,
  std.error_direct,
  p.value_bbs,
  p.value_fdr_bbs,
  p.value_fdr_direct,
  p.value_interaction_direct,
  p.value_interaction_fdr_direct,
  n_bbs,
  n_direct
)][, lapply(.SD, function(x) ifelse(is.na(x), "-", as.character(x)))],
rep_table, sep = "\t")

# GWAS IDs for Trial consistent metabolites — used to filter MR inputs
nmap <- read_xlsx(name_map_file, sheet = 1) |> as.data.table()
consistent_labels   <- replication[replicates %in% consistent_levels, .(label, platform)]
consistent_labels[nmap, gwas_id := i.gwas_id, on = c("label", "platform")]
dir.create(dirname(consistent_ids), recursive = TRUE, showWarnings = FALSE)
writeLines(na.omit(consistent_labels$gwas_id), consistent_ids)
tee("Consistent GWAS IDs written:", length(na.omit(consistent_labels$gwas_id)), "\n")
