#### DESCRIPTION #######################################
# Purpose of script:
# Take in the association statistics and plot volcano plots
# for the different studies
#
# Author: Nick Sunderland
#
# Date Created: 2025-04-17
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
assoc           <- snakemake@input[["assoc"]]
name_map        <- snakemake@input[["name_map"]]
bmi_mr          <- snakemake@input[["bmi_mr"]]
# associations
direct_dat      <- snakemake@output[["direct_dat"]]
bbs_dat         <- snakemake@output[["bbs_dat"]]
vol_bbs         <- snakemake@output[["vol_bbs"]]
vol_direct      <- snakemake@output[["vol_direct"]]
plat_corr       <- snakemake@output[["plat_corr"]]
plat_bland_alt  <- snakemake@output[["plat_bland_alt"]]
# trial & BMI-MR replications
corr_plot       <- snakemake@output[["corr_plot"]]
study_bland_alt <- snakemake@output[["study_bland_alt"]]
rep_plot        <- snakemake@output[["rep_plot"]]
rep_table       <- snakemake@output[["rep_table"]]
mr_rep_plot     <- snakemake@output[["mr_rep_plot"]]
# circos plot
circos_plot     <- snakemake@output[["circos_plot"]]
########################################################

# requirements
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(plotly)
library(htmlwidgets)
library(mcr)
library(readxl)
suppressPackageStartupMessages(library(data.table))


#############
# WEIGHT LOSS TRIAL ASSOCIATIONS
#############

# read
res <- lapply(assoc, fread) |> rbindlist(use.names=TRUE, fill=TRUE)


# grab number of metabolites
res[grepl("^olink_", feature_id), platform := "olink"]
total_metabs_assessed <- res[, .(n=uniqueN(label)), by=.(main_study, platform)]
rbind(
  total_metabs_assessed[, .(platform="all", n=sum(n)), by="main_study"],
  total_metabs_assessed
)


# subset to term of interest
res <- res[term == "timepointend" & model_name == "linear_mixed_time"]


# clean factors
res[, `:=`(main_study = as.factor(main_study),
           platform  = as.factor(platform),
           pathway_group = factor(pathway_group, levels = c("Lipid metabolism",
                                                            "Lipoproteins",
                                                            "Amino acid metabolism",
                                                            "Carbohydrate/energy metabolism",
                                                            "Nucleotide metabolism",
                                                            "Cofactors and vitamins",
                                                            "Inflammation",
                                                            "Fluid balance",
                                                            "Xenobiotics/other",
                                                            "Unknown")))]

# name_labels for top hits
n_top <- 10
res   <- res[order(p.value)]
res[, top_label := ifelse(seq_len(.N) < n_top, label, NA_character_), by=.(main_study, platform)]


# save direct and bbs associations
cols_to_save <- c("label", "pathway_group", "term", "estimate", "std.error", "statistic", "df", "p.value", "p.value_fdr", "p.value_holm", "p.value_interaction", "p.value_interaction_fdr", "p.value_interaction_holm",  "n", "data_type", "platform", "formula")
direct_assoc <- res[main_study == "direct", .SD, .SDcols = cols_to_save]
bbs_assoc    <- res[main_study == "bbs", .SD, .SDcols = cols_to_save]
fwrite(direct_assoc, direct_dat, sep="\t")
fwrite(bbs_assoc, bbs_dat, sep="\t")


# print out number of associations
cat(sprintf("DiRECT significant associations:
            total=%i; n=%i (%.1f%%) [nominal, p<0.05]; n=%i (%.1f%%) [fdr-ad, p<0.05]; n=%i (%.1f%%) [holm-ad, p<0.05],
            Nominal p<0.05: increasing=%i; decreasing=%i
            FDR p<0.05: increasing=%i; decreasing=%i
            HOLM p<0.05: increasing=%i; decreasing=%i \n",
            total_metabs_assessed[main_study=="direct", sum(n)],
            direct_assoc[!is.na(estimate), sum(p.value_interaction <0.05)],
            direct_assoc[!is.na(estimate), 100*sum(p.value_interaction<0.05) / total_metabs_assessed[main_study=="direct", sum(n)]],
            direct_assoc[!is.na(estimate), sum(p.value_interaction_fdr<0.05)],
            direct_assoc[!is.na(estimate), 100*sum(p.value_interaction_fdr<0.05) / total_metabs_assessed[main_study=="direct", sum(n)]],
            direct_assoc[!is.na(estimate), sum(p.value_interaction_holm<0.05)],
            direct_assoc[!is.na(estimate), 100*sum(p.value_interaction_holm<0.05) / total_metabs_assessed[main_study=="direct", sum(n)]],
            direct_assoc[!is.na(estimate) & p.value_interaction<0.05, sum(sign(estimate)==1)],
            direct_assoc[!is.na(estimate) & p.value_interaction<0.05, sum(sign(estimate)==-1)],
            direct_assoc[!is.na(estimate) & p.value_interaction_fdr<0.05, sum(sign(estimate)==1)],
            direct_assoc[!is.na(estimate) & p.value_interaction_fdr<0.05, sum(sign(estimate)==-1)],
            direct_assoc[!is.na(estimate) & p.value_interaction_holm<0.05, sum(sign(estimate)==1)],
            direct_assoc[!is.na(estimate) & p.value_interaction_holm<0.05, sum(sign(estimate)==-1)]))
cat(sprintf("BBS significant associations:
            total=%i; n=%i (%.1f%%) [nominal, p<0.05]; n=%i (%.1f%%) [fdr-ad, p<0.05]; n=%i (%.1f%%) [holm-ad, p<0.05],
            Nominal p<0.05: increasing=%i; decreasing=%i
            FDR p<0.05: increasing=%i; decreasing=%i
            HOLM p<0.05: increasing=%i; decreasing=%i \n",
            total_metabs_assessed[main_study=="bbs", sum(n)],
            bbs_assoc[!is.na(estimate), sum(p.value<0.05)],
            bbs_assoc[!is.na(estimate), 100*sum(p.value<0.05) / total_metabs_assessed[main_study=="bbs", sum(n)]],
            bbs_assoc[!is.na(estimate), sum(p.value_fdr<0.05)],
            bbs_assoc[!is.na(estimate), 100*sum(p.value_fdr<0.05) / total_metabs_assessed[main_study=="bbs", sum(n)]],
            bbs_assoc[!is.na(estimate), sum(p.value_holm<0.05)],
            bbs_assoc[!is.na(estimate), 100*sum(p.value_holm<0.05) / total_metabs_assessed[main_study=="bbs", sum(n)]],
            bbs_assoc[!is.na(estimate) & p.value<0.05, sum(sign(estimate)==1)],
            bbs_assoc[!is.na(estimate) & p.value<0.05, sum(sign(estimate)==-1)],
            bbs_assoc[!is.na(estimate) & p.value_fdr<0.05, sum(sign(estimate)==1)],
            bbs_assoc[!is.na(estimate) & p.value_fdr<0.05, sum(sign(estimate)==-1)],
            bbs_assoc[!is.na(estimate) & p.value_holm<0.05, sum(sign(estimate)==1)],
            bbs_assoc[!is.na(estimate) & p.value_holm<0.05, sum(sign(estimate)==-1)]))


# plot volcanos
volcanos <- list(bbs = vol_bbs, direct = vol_direct)
for (s in names(volcanos)) {
  df <- res[main_study == s, ]

  if (s == "bbs") {
    df$yval <- -log10(df$p.value)
  } else {
    df$yval <- -log10(df$p.value_interaction)
  }

  p <- ggplot(df, aes(x = estimate, y = yval)) +
    geom_point(aes(fill = pathway_group), shape = 21, color = "black", size = 2, stroke = 0.3, alpha = 0.7) +
    geom_text_repel(aes(label = top_label), size=2.5, color="black") +
    theme_light(base_size = 12) +
    labs(y = ifelse(s == "bbs", expression(-log[10](p.value)), expression(-log[10](p.value~interaction))),
         x = expression(beta[timepoint]),
         fill = "Pathway") +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(t = 30, r = 5, b = 5, l = 5)
    ) +
    facet_wrap(~platform,
               ncol = 1,
               scales = "free",
               labeller = labeller(platform = function(x) ifelse(x=="nightingale", "NMR", "MS")))

  # save
  ggsave(volcanos[[s]], p, width = 8, height = 8, dpi=300, bg="white")
}


#############
# CORRELATION BETWEEN PLATFORMS
#############
in_both_platforms <- res[, count_name := .N, by = c("label", "main_study")][count_name==2]
in_both_platforms <- dcast(in_both_platforms, label + main_study ~ platform, value.var = c("estimate", "std.error"))

# Deming regression line
deming_dat <- in_both_platforms[, {
  error_ratio <- mean(std.error_metabolon^2, na.rm=TRUE) / mean(std.error_nightingale^2, na.rm=TRUE)
  fit         <- mcr::mcreg(estimate_metabolon, estimate_nightingale, method.reg = "Deming", error.ratio = error_ratio)
  slope       <- coef(fit)[2]
  intercept   <- coef(fit)[1]
  cor_ci      <- confintr::ci_cor(estimate_metabolon, estimate_nightingale)
  pearson     <- cor_ci[["estimate"]]
  pearson_lb  <- cor_ci[["interval"]][[1]]
  pearson_ub  <- cor_ci[["interval"]][[2]]

  .(slope = coef(fit)[2], intercept = coef(fit)[1], pearsons = pearson, pearson_lb, pearson_ub,
    lab_x = -1,
    lab_y = ifelse(.BY[["main_study"]]=="bbs", 0.5, 0.4),
    label = sprintf("Pearson's R: %.2f (%.2f-%.2f) (%s)", pearson, pearson_lb, pearson_ub, toupper(.BY[["main_study"]])))
}, by="main_study"]

p <- ggplot(in_both_platforms, aes(x = estimate_metabolon, y = estimate_nightingale, color = main_study)) +
  geom_errorbar(aes(ymin = estimate_nightingale - 1.96*std.error_nightingale, ymax = estimate_nightingale + 1.96*std.error_nightingale), color = "lightgray", width=0) +
  geom_errorbarh(aes(xmin = estimate_metabolon - 1.96*std.error_metabolon, xmax = estimate_metabolon + 1.96*std.error_metabolon), color = "lightgray", height=0) +
  geom_abline(data = deming_dat, aes(intercept = intercept, slope = slope, color = main_study)) +
  geom_text_repel(aes(label = label), color="black") +
  geom_point() +
  geom_text(data = deming_dat, aes(x = lab_x, y = lab_y, label = label), fontface = "bold", hjust=0, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("bbs" = "#8DA0CB", "direct" = "#FC8D62"), labels = function(x) toupper(x)) +
  labs(x = "MS estimate", y = "NMR estimate", color = "Study")
p

# save
ggsave(plat_corr, p, width = 8, height = 6, dpi=300, bg="white")


# Bland-Altman
in_both_platforms[, `:=`(
  mean_estimate = (estimate_metabolon + estimate_nightingale) / 2,
  diff_estimate = estimate_metabolon - estimate_nightingale
)]
ba_summary_stats <- in_both_platforms[, .(` ` = 0,
                                          `Bias`        = mean(diff_estimate),
                                          `Upper limit of agreement` = mean(diff_estimate) + 1.96 * sd(diff_estimate),
                                          `Lower limit of agreement` = mean(diff_estimate) - 1.96 * sd(diff_estimate))]
ba_summary_stats <- melt(ba_summary_stats)
ba_summary_stats[, `:=`(linetype = ifelse(variable == " ", "solid", "dashed"))]
ba_summary_stats

# report number outside 95% CI
in_both_platforms[diff_estimate > ba_summary_stats[variable=="Upper limit of agreement", value] |
                  diff_estimate < ba_summary_stats[variable=="Lower limit of agreement", value],
                  .(num_outlier = .N,
                    pct_outlier = .N / nrow(in_both_platforms))]

p <- ggplot(in_both_platforms, aes(x = mean_estimate, y = diff_estimate, fill = main_study)) +
  geom_hline(data = ba_summary_stats, aes(yintercept = value, linetype = linetype)) +
  geom_point(shape = 21, col = "black", size = 3) +
  geom_text(data = ba_summary_stats,
            aes(label = variable, x = max(in_both_platforms$mean_estimate, na.rm=T)+.02, y = value+.02),
            hjust=1,
            color="gray",
            inherit.aes = F,
            size = 3) +
  geom_text_repel(aes(label = label), color = "black", size = 3) +
  scale_fill_manual(values = c("bbs" = "#8DA0CB", "direct" = "#FC8D62"), labels = function(x) toupper(x)) +
  scale_linetype_identity() +
  theme_light() +
  theme(legend.position = "top") +
  labs(y = "Difference",
       x = "Mean",
       fill = "Study")
p

# save
ggsave(plat_bland_alt, p, width = 8, height = 6, dpi=300, bg="white")


#############
# REPLICATION
#############

# pivot wider to assess replication
replication <- dcast(res, label + platform ~ main_study, value.var = c("estimate", "std.error", "p.value", "p.value_fdr", "p.value_interaction", "p.value_interaction_fdr", "p.value_holm", "n"))

replication[, `:=`(coverage         = fcase(!is.na(estimate_direct) & !is.na(estimate_bbs), "both",
                                            is.na(estimate_direct) & !is.na(estimate_bbs), "bbs_only",
                                            !is.na(estimate_direct) & is.na(estimate_bbs), "direct_only",
                                            default = "none"),
                   same_direction   = sign(estimate_bbs) == sign(estimate_direct),
                   both_significant = p.value_fdr_bbs<0.05 & p.value_interaction_fdr_direct<0.05)]


replication[, replicates := factor(fcase(same_direction & both_significant, "Trial consistent",
                                         both_significant,                  "Trial FDR P < 0.05 only",
                                         same_direction,                    "Trial direction only",
                                         default = "None"),
                                   levels = c("Trial consistent", "Trial FDR P < 0.05 only", "Trial direction only", "None"))]

setkeyv(replication, c("replicates", "label"))

# print overlap
non_overlapping <- replication[is.na(estimate_bbs) | is.na(estimate_direct), ]

replication[!is.na(estimate_bbs) & !is.na(estimate_direct), .N, by=platform]

cat(sprintf("Total features=%i (BBS=%i; DiRECT=%i), Overlap n=%i (%.1f%%), concordant=%i (%.1f%%), disconcordant=%i (%.1f%%)",
            replication[, uniqueN(interaction(label, platform))],
            replication[!is.na(estimate_bbs), uniqueN(interaction(label, platform))],
            replication[!is.na(estimate_direct), uniqueN(interaction(label, platform))],
            replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
            100 * replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)] / replication[, uniqueN(label)],
            replication[, sum(same_direction, na.rm = TRUE)],
            100 * replication[, sum(same_direction, na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
            replication[, sum(!same_direction, na.rm = TRUE)],
            100 * replication[, sum(!same_direction, na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)]))

cat(sprintf("Replication: sign=%i (%.1f%%), sig=%i (%.1f%%), both=%i (%.1f%%)
            both-positive direction=%i (%.1f%%); both-negative direction=%i (%.1f%%)",
            replication[, sum(same_direction, na.rm = TRUE)],
            100 * replication[, sum(same_direction, na.rm = TRUE)] / replication[, sum(!is.na(same_direction))],
            replication[, sum(both_significant, na.rm = TRUE)],
            100 * replication[, sum(both_significant, na.rm = TRUE)] / replication[, sum(!is.na(both_significant))],
            replication[, sum(replicates=="Trial consistent", na.rm = TRUE)],
            100 * replication[, sum(replicates=="Trial consistent", na.rm = TRUE)] / replication[, sum(!is.na(estimate_bbs) & !is.na(estimate_direct), na.rm = TRUE)],
            replication[, sum(replicates=="Trial consistent" & sign(estimate_bbs)==1, na.rm = TRUE)],
            100 * replication[, sum(replicates=="Trial consistent" & sign(estimate_bbs)==1, na.rm = TRUE)] / replication[, sum(replicates=="Trial consistent", na.rm = TRUE)],
            replication[, sum(replicates=="Trial consistent" & sign(estimate_bbs)==-1, na.rm = TRUE)],
            100 * replication[, sum(replicates=="Trial consistent" & sign(estimate_bbs)==-1, na.rm = TRUE)] / replication[, sum(replicates=="Trial consistent", na.rm = TRUE)]))

top_changers <- replication[replicates == "Trial consistent",
                            .SD[unique(c(
                              head(order(-abs(estimate_bbs)), 5),
                              head(order(-abs(estimate_direct)), 5)
                            ))]
                            ]
top_changers



#############
# CORRELATION
#############
corr_dat    <- dcast(res[model_name=="linear_mixed_time"][abs(estimate) > 2, estimate := NA_real_],
                     label + platform ~ main_study, value.var = c("estimate", "std.error"))
corr_dat    <- corr_dat[complete.cases(corr_dat)]

# Deming regression line
deming_dat <- corr_dat[, {
  error_ratio <- mean(std.error_bbs^2, na.rm=TRUE) / mean(std.error_direct^2, na.rm=TRUE)
  fit         <- mcr::mcreg(estimate_bbs, estimate_direct, method.reg = "Deming", error.ratio = error_ratio)
  slope       <- coef(fit)[2]
  intercept   <- coef(fit)[1]
  cor_ci      <- confintr::ci_cor(estimate_bbs, estimate_direct)
  pearson     <- cor_ci[["estimate"]]
  pearson_lb  <- cor_ci[["interval"]][[1]]
  pearson_ub  <- cor_ci[["interval"]][[2]]

  .(slope = coef(fit)[2], intercept = coef(fit)[1], pearson, pearson_lb, pearson_ub,
    lab_x = -1.2,
    lab_y = 1.15,
    label = sprintf("Pearson's R: %.2f (%.2f-%.2f)", pearson, pearson_lb, pearson_ub))
}, by="platform"]

# label greatest orthogonal residuals
corr_dat[deming_dat, `:=`(slope=i.slope, intercept=i.intercept), on="platform"]
corr_dat[, orth_resid := abs(slope * estimate_bbs - estimate_direct + intercept) / sqrt(slope^2 + 1)]
corr_dat[, labels_to_plot := fifelse(frank(-orth_resid, ties.method = "first") <= 10, label, NA_character_), by="platform"]

p <- ggplot(corr_dat, aes(x = estimate_bbs, y = estimate_direct, color = platform)) +
  geom_point() +
  geom_abline(data = deming_dat, aes(intercept = intercept, slope = slope, color = platform)) +
  geom_text_repel(aes(label = labels_to_plot), color="black", size = 2.5, na.rm = TRUE) +
  geom_text(data = deming_dat, aes(x = lab_x, y = lab_y, label = label), fontface = "bold", hjust=0, show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  theme_light(base_size = 14) +
  labs(x = expression(beta[timepoint~BBS]), y = expression(beta[timepoint~DiRECT])) +
  guides(color = "none") +
  facet_wrap(~platform, ncol=1, labeller = labeller(platform = function(x) ifelse(x=="nightingale", "NMR", "MS")))
p

# save
ggsave(corr_plot, p, width = 7, height = 10, dpi=300)



# Bland-Altman
corr_dat[, `:=`(
  mean_estimate = (estimate_bbs + estimate_direct) / 2,
  diff_estimate = estimate_bbs - estimate_direct
)]
ba_study_summary_stats <- corr_dat[, .(` ` = 0,
                                       `Bias`        = mean(diff_estimate),
                                       `Upper limit of agreement` = mean(diff_estimate) + 1.96 * sd(diff_estimate),
                                       `Lower limit of agreement` = mean(diff_estimate) - 1.96 * sd(diff_estimate))]
ba_study_summary_stats <- melt(ba_study_summary_stats)
ba_study_summary_stats[, `:=`(linetype = ifelse(variable == " ", "solid", "dashed"))]
ba_study_summary_stats

# report number outside 95% CI
corr_dat[diff_estimate > ba_study_summary_stats[variable=="Upper limit of agreement", value] |
           diff_estimate < ba_study_summary_stats[variable=="Lower limit of agreement", value], .(num_outlier = .N,
                                                                                                  pct_outlier = .N / nrow(corr_dat))]

p <- ggplot(corr_dat, aes(x = mean_estimate, y = diff_estimate, fill = platform)) +
  geom_hline(data = ba_study_summary_stats, aes(yintercept = value, linetype = linetype)) +
  geom_point(shape = 21, stroke=0, size = 2, alpha=0.7) +
  geom_text(data = ba_study_summary_stats,
            aes(label = variable, x = max(corr_dat$mean_estimate, na.rm=T), y = value),
            hjust=1, vjust=-0.5,
            color="gray",
            inherit.aes = F,
            size = 3) +
  geom_text_repel(data = corr_dat[diff_estimate > ba_study_summary_stats[variable=="Upper limit of agreement", value] |
                                    diff_estimate < ba_study_summary_stats[variable=="Lower limit of agreement", value]],
                  aes(label = label), color = "black", size = 3) +
  scale_fill_brewer(palette = "Set1", labels = function(x) ifelse(x=="nightingale", "NMR", "MS")) +
  scale_linetype_identity() +
  theme_light() +
  theme(legend.position = "top") +
  labs(y = "Difference",
       x = "Mean",
       fill = "Platform")
p

# save
ggsave(study_bland_alt, p, width = 10, height = 7, dpi=300, bg="white")



# plot replication
deming_dat <- replication[, {
  error_ratio <- mean(std.error_bbs^2, na.rm=TRUE) / mean(std.error_direct^2, na.rm=TRUE)
  fit         <- mcr::mcreg(estimate_bbs, estimate_direct, method.reg = "Deming", error.ratio = error_ratio, na.rm = TRUE)
  slope       <- coef(fit)[2]
  intercept   <- coef(fit)[1]
  .(slope = coef(fit)[2], intercept = coef(fit)[1])
}]

p <- ggplot(replication,
            aes(x = estimate_bbs, y = estimate_direct, color = replicates, shape = replicates)) +
  geom_abline(data = deming_dat, aes(intercept = intercept, slope = slope), color = "black") +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "darkgray", linetype="dotted") +
  geom_vline(xintercept = 0, color = "darkgray", linetype="dotted") +
  geom_text_repel(data = replication[replicates=="Trial consistent" | replicates=="Trial FDR P < 0.05 only"],
                  aes(label = label), size=3, color="black") +
  scale_color_manual(values = c("Trial consistent"="red", "Trial FDR P < 0.05 only"=scales::muted("blue"), "Trial direction only"="gray", "None"="lightgray")) +
  theme_classic(base_size = 14) +
  labs(x = expression(beta[timepoint~BBS]), y = expression(beta[timepoint~DiRECT])) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(t = 30, r = 5, b = 5, l = 5)
  )
p

# save
ggsave(rep_plot, p, width = 10, height = 10, dpi=300)



#######################
# BMI adjustment of replicating metabolites
#######################
bbs_bmi_change <- -10.5
bbs_bmi_sd <- 7.3

replication[replicates=="Trial consistent" ,
            `:=`(bmi_estimate_bbs    = estimate_bbs / (bbs_bmi_change/bbs_bmi_sd),
                 bmi_std.error_bbs   = std.error_bbs / abs(bbs_bmi_change/bbs_bmi_sd),
                 bmi_p.value_bbs     = p.value_bbs,
                 bmi_p.value_fdr_bbs = p.value_fdr_bbs)]

# -4 for intervention, -1 for control --> 0
direct_bmi_change <- -3.7
direct_bmi_sd <- 4.5
replication[replicates=="Trial consistent" ,
            `:=`(bmi_estimate_direct    = estimate_direct / (direct_bmi_change/direct_bmi_sd),
                 bmi_std.error_direct   = std.error_direct / abs(direct_bmi_change/direct_bmi_sd),
                 bmi_p.value_direct     = p.value_direct,
                 bmi_p.value_fdr_direct = p.value_fdr_direct,
                 bmi_p.value_interaction_direct = p.value_interaction_direct,
                 bmi_p.value_interaction_fdr_direct = p.value_interaction_fdr_direct)]



#######################
# BMI metabolite MR results
#######################
map <- read_xlsx(name_map, sheet=1) |> as.data.table()
mr <- fread(bmi_mr)
mr[map, `:=`(label = i.label), on=c("outcome"="gwas_name", "metabolite_platform"="platform")]

# take IVW MR result from most significant platform / GWAS
mr <- mr[method=="mr_ivw", .SD[which.min(p)], by = c("label", "metabolite_platform")]
#mr[, replicates := "Trial consistent"]

# add mr results
replication[mr[method=="mr_ivw" & label!=""],
            `:=`(mr_nsnp             = i.n_snp,
                 metabolite_gwas_platform = i.metabolite_platform,
                 metabolite_gwas_source = fcase(i.metabolite_platform=="nightingale", "karjalainen_2024",
                                                i.metabolite_platform=="metabolon", "chen_2023", default = NA_character_),
                 bmi_estimate_mr     = i.b,
                 bmi_std.error_mr    = i.b_se,
                 bmi_p.value_mr      = i.p,
                 bmi_p.value_fdr_mr  = p.adjust(i.p, method="fdr")),
             on=c("label")]
replication[map, metabolite_gwas_id := i.gwas_id, on = c("label", "metabolite_gwas_platform"="platform")]
replication[, metabolite_gwas_platform := NULL]

# annotate replication
replication[, `:=`(all_same_direction = sign(bmi_estimate_bbs) == sign(bmi_estimate_direct) & sign(bmi_estimate_bbs) == sign(bmi_estimate_mr),
                   all_significant    = bmi_p.value_fdr_bbs < 0.05 & bmi_p.value_interaction_fdr_direct < 0.05 & bmi_p.value_fdr_mr < 0.05)]

replication[, replicates_mr := factor(fcase(all_same_direction & all_significant,                       "Trial & MR consistent",
                                            replicates=="Trial consistent" & bmi_p.value_fdr_mr < 0.05, "Trial consistent & MR FDR P < 0.05 only",
                                            replicates=="Trial consistent" & all_same_direction,        "Trial consistent & MR direction only",
                                            !is.na(bmi_estimate_mr), "None",
                                            default = NA_character_),
                                      levels = c("Trial & MR consistent", "Trial consistent & MR FDR P < 0.05 only", "Trial consistent & MR direction only", "None"))]



# print replication results
replication[map, pathway := i.pathway_group, on=c("label","platform")]
replication[, label_platform := interaction(label, platform)]

replication[replicates_mr=="Trial & MR consistent", .(num = uniqueN(label_platform), pct = uniqueN(label_platform) / replication[replicates_mr=="Trial & MR consistent", uniqueN(label_platform)]), by=c("pathway")][order(-pct)]

setcolorder(replication, "label_platform")
rbind(
  replication[, .(cat = "Total metabolites", num = uniqueN(label_platform))],
  replication[!is.na(estimate_bbs), .(cat = "Total BBS", num = uniqueN(label_platform))],
  replication[, .(cat = "BBS num FDR[all]", num = sum(p.value_fdr_bbs<0.05, na.rm=TRUE))],
  replication[, .(cat = "BBS num FDR[+ve with WL]", num = sum(p.value_fdr_bbs<0.05 & estimate_bbs>0, na.rm=TRUE))],
  replication[, .(cat = "BBS num FDR[-ve with WL]", num = sum(p.value_fdr_bbs<0.05 & estimate_bbs<0, na.rm=TRUE))],
  replication[!is.na(estimate_direct), .(cat = "Total DIR", num = uniqueN(label_platform))],
  replication[, .(cat = "DIR num FDR[all]", num = sum(p.value_interaction_fdr_direct<0.05, na.rm=TRUE))],
  replication[, .(cat = "DIR num FDR[+ve with WL]", num = sum(p.value_interaction_fdr_direct<0.05 & estimate_direct>0, na.rm=TRUE))],
  replication[, .(cat = "DIR num FDR[-ve with WL]", num = sum(p.value_interaction_fdr_direct<0.05 & estimate_direct<0, na.rm=TRUE))],
  replication[, .(cat = "Total in BBS & DIR", num = sum(!is.na(estimate_direct) & !is.na(estimate_bbs)))],
  setnames(replication[, .(num = uniqueN(label_platform)), by="replicates"], "replicates", "cat"),
  replication[, .(cat = "Trial consistent[+ve with WL]", num = sum(replicates=="Trial consistent" & estimate_direct>0))],
  replication[, .(cat = "Trial consistent[-ve with WL]", num = sum(replicates=="Trial consistent" & estimate_direct<0))],
  replication[, .(cat = "Trial consistent & any MR result", num = sum(replicates=="Trial consistent" & !is.na(bmi_estimate_mr)))],
  replication[, .(cat = "Trial & MR FDR<0.05", num = sum(replicates=="Trial consistent" & bmi_p.value_fdr_mr<0.05, na.rm=T))],
  replication[, .(cat = "Trial & MR dir consistent & FDR<0.05", num = sum(replicates=="Trial consistent" & sign(bmi_estimate_mr)==sign(bmi_estimate_bbs) & bmi_p.value_fdr_mr<0.05, na.rm=T))],
  replication[, .(cat = "Trial & MR dir inconsistent & FDR<0.05", num = sum(replicates=="Trial consistent" & sign(bmi_estimate_mr)!=sign(bmi_estimate_bbs) & bmi_p.value_fdr_mr<0.05, na.rm=T))],
  setnames(replication[, .(num = uniqueN(label_platform)), by="replicates_mr"], "replicates_mr", "cat")
)

# save results
fwrite(replication[, .(
  label,
  platform,
  replicates_trials       = replicates,
  replicates_trials_mr    = replicates_mr,
  time_estimate_bbs       = estimate_bbs,
  time_estimate_direct    = estimate_direct,
  time_std.error_bbs      = std.error_bbs,
  time_std.error_direct   = std.error_direct,
  time_p.value_bbs        = p.value_bbs,
  time_p.value_direct     = p.value_direct,
  time_p.value_fdr_bbs    = p.value_fdr_bbs,
  time_p.value_fdr_direct = p.value_fdr_direct,
  time_p.value_interaction_fdr_direct = p.value_interaction_fdr_direct,
  n_bbs                   = n_bbs,
  n_direct                = n_direct,
  metabolite_gwas_id, mr_nsnp,
  bmi_estimate_bbs, bmi_estimate_direct, bmi_estimate_mr,
  bmi_std.error_bbs, bmi_std.error_direct, bmi_std.error_mr,
  bmi_p.value_bbs, bmi_p.value_direct, bmi_p.value_interaction_direct,  bmi_p.value_mr,
  bmi_p.value_fdr_bbs, bmi_p.value_fdr_direct, bmi_p.value_interaction_fdr_direct, bmi_p.value_fdr_mr
)][, lapply(.SD, function(x) ifelse(is.na(x), "-", as.character(x)))], rep_table, sep = "\t")




# CIRCOS PLOT
library(circlize)

# --- 1. Data ---
df <- replication

# --- 2. Assign Colors ---
alpha <- 0.7
alpha2 <- 0.3
color_list <- list(
  "Consistent -ve BMI association" = "firebrick",                                # 1 - decreases with BMI
  "Consistent +ve BMI association" = "royalblue",                                # 2 - increases with BMI
  "-ve BMI association"            = adjustcolor("firebrick", alpha.f = alpha2), # 3 - decreases with BMI / increases with WL
  "+ve BMI association"            = adjustcolor("royalblue", alpha.f = alpha2), # 4 - increases with BMI / decreases with WL
  "No BMI association"             = adjustcolor("lightgrey", alpha.f = alpha2), # 5 - no associated
  "No data"                        = "white"                                     # 6 - no result
)
df[, bbs_color := fcase(as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_bbs < 0, 1,
                        as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_bbs > 0, 2,
                        p.value_fdr_bbs < 0.05 & estimate_bbs < 0,                                     4,
                        p.value_fdr_bbs < 0.05 & estimate_bbs > 0,                                     3,
                        is.na(estimate_bbs),                                                           6,
                        default = 5)]
df[, dir_color := fcase(as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_direct < 0, 1,
                        as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_direct > 0, 2,
                        p.value_fdr_direct < 0.05 & estimate_direct < 0,                                         4,
                        p.value_fdr_direct < 0.05 & estimate_direct > 0,                                         3,
                        is.na(estimate_direct),                                                           6,
                        default = 5)]
df[, mr_color  := fcase(as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_mr < 0, 1,
                        as.character(replicates_mr) == "Trial & MR consistent" & bmi_estimate_mr > 0, 2,
                        bmi_p.value_fdr_mr < 0.05 & bmi_estimate_mr < 0,                              3,
                        bmi_p.value_fdr_mr < 0.05 & bmi_estimate_mr > 0,                              4,
                        is.na(bmi_estimate_mr),                                                       6,
                        default = 5)]
df[, `:=`(
  bbs_color = unlist(color_list[bbs_color]),
  dir_color = unlist(color_list[dir_color]),
  mr_color  = unlist(color_list[mr_color]),
  color_order = factor(unlist(color_list[dir_color]), levels = unname(color_list))
)]

df <- df[order(color_order)]


# --- 4. Initialize Circos ---
{
  png(filename = circos_plot,
      width = 2000, height = 1600, res = 300)
  par(xpd = TRUE, mar = c(1, 1, 1, 1))

  circos.clear()
  circos.par(start.degree = -5,
             gap.degree   = 0,
             cell.padding = c(0, 0, 0, 0),
             canvas.xlim = c(-1.0, 2.2),
             canvas.ylim = c(-1.2, 1.2)
             )

  n_metabolites <- nrow(df)
  circos.initialize(factors = df$label_platform,
                    xlim = cbind(rep(0, n_metabolites), rep(1, n_metabolites)))

  # --- 5. Add Rings for Each Trial ---
  n_rings <- 3
  label_indices <- which(df$dir_color %in% c("firebrick", "royalblue"))


  # Trial 3 ring
  circos.trackPlotRegion(track.index = n_rings-2,
                         ylim = c(0, 1),
                         bg.border = NA,
                         panel.fun = function(x, y) {
                           i <- CELL_META$sector.numeric.index
                           col <- df$mr_color[i]
                           circos.rect(0, 0, 1, 1, col = col, border = NA)
                         })


  # Trial 2 ring
  circos.trackPlotRegion(track.index = n_rings-1, ylim = c(0, 1), bg.border = NA,
                         panel.fun = function(x, y) {
                           i <- CELL_META$sector.numeric.index
                           col <- df$bbs_color[i]
                           circos.rect(0, 0, 1, 1, col = col, border = NA)
                         })

  # Trial 1 ring
  circos.trackPlotRegion(track.index = n_rings-0, ylim = c(0, 1), bg.border = NA,
                         panel.fun = function(x, y) {
                           i <- CELL_META$sector.numeric.index
                           col <- df$dir_color[i]
                           circos.rect(0, 0, 1, 1, col = col, border = NA)
                         })

  legend("topleft",
         legend = names(color_list),
         fill   = unname(unlist(color_list)),
         border = "black",
         bty    = "n",
         cex    = 0.9,
         y.intersp = 1.2,
         x.intersp = 0.3,
         pt.cex      = 1.9,
         seg.len     = 1.9,
         inset  = c(0.62, 0.12))

  dev.off()
}


# subset to the trial replicating metabolites
trial_and_mr_replicating <- melt(replication[replicates=="Trial consistent"],
                                 id.vars = c("label", "platform", "replicates_mr"),
                                 variable.name = "method",
                                 measure.vars = list(
                                   estimate = grep("^bmi_estimate_(bbs|direct|mr)$", names(replication), value = TRUE),
                                   se       = grep("^bmi_std\\.error_(bbs|direct|mr)$", names(replication), value = TRUE),
                                   p        = grep("^bmi_p\\.value_(bbs|direct|mr)$", names(replication), value = TRUE),
                                   p_fdr    = grep("^bmi_p\\.value_fdr_(bbs|direct|mr)$", names(replication), value = TRUE)
                                 ))

trial_and_mr_replicating[, method := factor(method, levels = 1:3, labels = c("bbs", "direct", "mr"))]
trial_and_mr_replicating[, type := factor(method, levels = c("bbs", "direct", "mr"), labels = c("wl", "wl", "mr"))]
trial_and_mr_replicating[, label := factor(label, levels = unique(label[order(-estimate)]))]

trial_and_mr_replicating <- dcast(trial_and_mr_replicating[!is.na(method)], label + platform + replicates_mr + method ~ type, value.var = c("estimate", "se", "p", "p_fdr"))
trial_and_mr_replicating[, `:=`(estimate_mr = stats::na.omit(estimate_mr),
                                se_mr       = stats::na.omit(se_mr),
                                p_mr        = stats::na.omit(p_mr),
                                p_fdr_mr    = stats::na.omit(p_fdr_mr)), by=c("label", "platform")]
trial_and_mr_replicating <- trial_and_mr_replicating[method != "mr" & !is.na(estimate_mr)]

# recheck number of metabs
trial_and_mr_replicating[, uniqueN(interaction(label, platform)), by="replicates_mr"]


# plot MR replication
p <- ggplot(trial_and_mr_replicating[abs(estimate_wl) < 2], aes(x = estimate_wl, y = estimate_mr, color = replicates_mr, shape = method)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0, color = "darkgray", linetype="dotted") +
  geom_vline(xintercept = 0, color = "darkgray", linetype="dotted") +
  geom_text_repel(data = trial_and_mr_replicating[abs(estimate_wl) < 2 & grepl("(Trial & MR consistent)|(Trial consistent & MR FDR P < 0.05 only)", replicates_mr)],
                  aes(label = label), size=3, color="black") +
  scale_color_manual(values = c("Trial & MR consistent"="red", "Trial consistent & MR FDR P < 0.05 only"=scales::muted("blue"), "Trial consistent & MR direction only"="gray", "None"="lightgray")) +
  scale_shape_discrete(labels = function(x) toupper(x)) +
  theme_classic(base_size = 14) +
  labs(y = "MR - metabolite change per +1 SD genetically predicted BMI change",
       x = "Metabolite change per +1 SD BMI change with weight loss intervention") +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(t = 30, r = 5, b = 5, l = 5)
  )
p

# save
ggsave(mr_rep_plot, p, width = 10, height = 10, dpi=300)
