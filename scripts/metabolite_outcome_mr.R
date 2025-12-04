#### DESCRIPTION #######################################
# Purpose of script:
# Take ...
#
# Author: Nick Sunderland
#
# Date Created: 2025-04-17
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
name_map                    <- snakemake@input[["name_map"]]
outcome_mr                  <- snakemake@input[["outcome_mr"]]
metab_loci_fp               <- snakemake@input[["metab_loci_fp"]]
replicating                 <- snakemake@input[["replicating"]]
outcome_snp_overlap_dir     <- snakemake@params[["outcome_snp_overlap_dir"]]
outcome_mr_tbl              <- snakemake@output[["outcome_mr_tbl"]]
outcome_mr_plot             <- snakemake@output[["outcome_mr_plot"]]
outcome_mr_forest           <- snakemake@output[["outcome_mr_forest"]]
outcome_volcano_forest      <- snakemake@output[["outcome_volcano_forest"]]
outcome_mr_sig              <- snakemake@output[["outcome_mr_sig"]]
outcome_instr_tbl           <- snakemake@output[["outcome_instr_tbl"]]
pathway_overlap_tbl         <- snakemake@output[["pathway_overlap_tbl"]]
bmi_metabolite_snps         <- snakemake@output[["bmi_metabolite_snps"]]
bmi_metabolite_variant_venn <- snakemake@output[["bmi_metabolite_variant_venn"]]
liu_model1_overlap_venn_ms  <- snakemake@output[["liu_model1_overlap_venn_ms"]]
liu_model1_overlap_venn_nmr <- snakemake@output[["liu_model1_overlap_venn_nmr"]]
out_mr_obs_heatmap          <- snakemake@output[["out_mr_obs_heatmap"]]
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


# read
outcomes <- fread(outcome_mr)
map      <- read_xlsx(name_map, sheet=1) |> as.data.table()
rep      <- fread(replicating)
rep[, label_plat := paste0(label, " [", platform, "]")]

# clean
outcomes <- outcomes[outcome %in% c("heart_failure", "hfpef", "hfref")]
outcomes[map, `:=`(label = i.label, pathway_group=i.pathway_group, gwas_id=i.gwas_id), on=c("exposure"="gwas_name")]
outcomes <- rep[replicates_trials_mr=="Trial & MR consistent", ][outcomes,  on=c("label"="label", "platform"="metabolite_platform"), nomatch=NULL]
outcomes[, `:=`(pathway_group = factor(pathway_group, levels = c("Lipid metabolism",
                                                                 "Lipoproteins",
                                                                 "Amino acid metabolism",
                                                                 "Carbohydrate/energy metabolism",
                                                                 "Nucleotide metabolism",
                                                                 "Cofactors and vitamins",
                                                                 "Inflammation",
                                                                 "Fluid balance",
                                                                 "Xenobiotics/other",
                                                                 "Unknown")))]
outcomes[, mean_nsnp := round(mean(n_snp, na.rm=T)), by=c("exposure", "platform")]
outcomes[, exposure_lab := paste0(exposure, " ", ifelse(platform=="nightingale","(N; ","(M; "), "nsnp=", mean_nsnp, ")")]
outcomes[, exposure_lab := factor(exposure_lab, levels = unique(exposure_lab[order(pathway_group)]))]
outcomes[, outcome := factor(outcome, levels = c("heart_failure","hfpef","hfref"), labels = c("All-cause HF", "Non-ischaemic HFpEF", "Non-ischaemic HFrEF"))]
outcomes <- outcomes[!is.na(outcome)]
outcomes[, effect_dir := ifelse(sign(b)==-1, "Neg", "Pos")]
outcomes[, method := factor(method, levels = rev(c("mr_ivw", "mr_egger", "mr_weighted_median", "mr_weighted_mode")))]

# save
fwrite(outcomes[, .(
  exposure = label,
  platform = platform,
  gwas_id  = gwas_id,
  pathway  = pathway_group,
  n_snp    = n_snp,
  fstat    = fstat,
  outcome  = outcome,
  mr_method= method,
  or_ci    = ifelse(!is.na(b), sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b-1.96*b_se), exp(b+1.96*b_se)), NA_character_),
  beta     = b,
  std_error= b_se,
  p_value  = p,
  intercept= intercept,
  int_se   = int_se,
  int_p    = int_p,
  qstat    = qstat,
  qstat_p  = qstat_p,
  `.`      = NA_character_,
  bbs_bmi_estimate       = bmi_estimate_bbs,
  direct_bmi_estimate    = bmi_estimate_direct,
  mr_bmi_estimate        = bmi_estimate_mr,
  bbs_bmi_std_error      = bmi_std.error_bbs,
  direct_bmi_std_error   = bmi_std.error_direct,
  mr_bmi_std_error       = bmi_std.error_mr,
  bbs_bmi_p_value        = bmi_p.value_bbs,
  direct_bmi_p_value     = bmi_p.value_direct,
  direct_bmi_p_value_interaction = bmi_p.value_interaction_direct,
  mr_bmi_p_value         = bmi_p.value_mr,
  bbs_bmi_fdr_p_value    = bmi_p.value_fdr_bbs,
  direct_bmi_fdr_p_value = bmi_p.value_fdr_direct,
  direct_bmi_interaction_fdr_p_value = bmi_p.value_interaction_fdr_direct,
  mr_bmi_fdr_p_value     = bmi_p.value_fdr_mr
)][, lapply(.SD, function(x) ifelse(x==""|is.na(x), ".", as.character(x)))], outcome_mr_tbl, sep="\t")


# report successful MR
any_mr_success <- outcomes[, .(any_mr_success = any(!is.na(n_snp))), by = interaction(label, platform)]
any_mr_success[, n_tot := .N]
any_mr_success[, .(n_tot = n_tot[1],
                   successful_instrument = .N,
                   successful_instrument_pct = .N / n_tot[1],
                   min_nsnp = min(outcomes$n_snp, na.rm=T),
                   max_nsnp = max(outcomes$n_snp, na.rm=T)), by = "any_mr_success"]
outcomes[, .(mean_snps = mean(n_snp, na.rm=T)), by=c("pathway_group")][order(-mean_snps)]


# plot
p <- ggplot(outcomes[method=="mr_ivw" & !grepl("^X", as.character(exposure))],
            aes(y = exposure_lab, x = outcome, fill = pathway_group, color=effect_dir, alpha = ifelse(p<0.05, -log10(p), 0))) +
  geom_tile(color=NA) +
  geom_point(data = outcomes[p<0.05 & method=="mr_ivw" & !grepl("^X", as.character(exposure)), ], alpha=1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  scale_y_discrete(labels = function(x) sub("(?i) levels", "", x)) +
  labs(
    color = "d"
  ) +
  theme(
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    alpha = "-log10(p)",
    fill = "Pathway",
    color = "P<0.05 / Effect direction",
    y = "Weight loss & BMI MR significant metabolite (exposure)",
    x = "Outcome"
  )
p

# save
ggsave(outcome_mr_plot, p, width = 6, height = 10, dpi=300, bg="white")



# Forest plot for HFpEF vs HFrEF
forest_dat <- outcomes[method=="mr_ivw"]
forest_dat[, `:=`(or = exp(b), or_lb = exp(b-1.96*b_se), or_ub = exp(b+1.96*b_se))]
forest_dat[, p_fdr := p.adjust(p, method="fdr", n=.N), by="outcome"]
forest_dat[, p_sig := factor(fcase(p_fdr < 0.05, 1,
                                   p     < 0.05, 2,
                                   default = 3), levels=1:3, labels=c("FDR P < 0.05", "P < 0.05", "P \u2265 0.05"))]
forest_dat[
  , exposure_lab := factor(
    exposure_lab,
    levels = rev(exposure_lab[grepl("(?i)all", outcome)][
      order(pathway_group[grepl("(?i)all", outcome)], -b[grepl("(?i)all", outcome)])
    ])
  )
]
p <- ggplot(forest_dat[!grepl("^X", as.character(exposure_lab))], aes(y = exposure_lab, x = or, color=pathway_group, )) +
  geom_vline(xintercept = 1, color="darkgray") +
  geom_errorbarh(aes(xmin = or_lb, xmax = ifelse(or_ub>2, 2, or_ub)), color="lightgray", height=0) +
  geom_point(aes(shape = p_sig)) +
  scale_color_brewer(palette = "Set1", direction = 1) +
  scale_shape_manual(values = c("P \u2265 0.05"=21, "P < 0.05"=19, "FDR P < 0.05" = 17)) +
  labs(x = "OR for outcome per 1 SD increase in metabolite",
       shape = "P-value",
       size = "P-value",
       color = "Pathway") +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size=6),
    axis.title.y = element_blank()
  ) +
  facet_wrap(~outcome, scales = "free_x")
p

# save
ggsave(outcome_mr_forest, p, width = 10, height = 12, dpi=300, bg="white")


# Volcano of all-cause HF results
library(cowplot)
mr_volcano_dat <- fread(outcome_mr)
mr_volcano_dat[map, `:=`(label = i.label, pathway_group=i.pathway_group, gwas_id=i.gwas_id), on=c("exposure"="gwas_name")]
mr_volcano_dat[, bmi_metabolite := factor("No", levels=c("No","Yes"))]
mr_volcano_dat[rep[replicates_trials_mr=="Trial & MR consistent"], bmi_metabolite := factor("Yes", levels=c("No","Yes")), on=c("label"="label", "metabolite_platform"="platform")]
mr_volcano_dat <- mr_volcano_dat[method=="mr_ivw" & as.character(outcome)=="heart_failure"]

p <- ggplot(mr_volcano_dat, aes(x = b, y = -log10(p), alpha=bmi_metabolite)) +
  geom_point(aes(color=bmi_metabolite), shape=19, size = 2) +
  geom_hline(yintercept = -log10(0.05), color="gray", linetype="dashed") +
  scale_color_manual(values = c("Yes" = "red", "No" = "royalblue")) +
  scale_alpha_manual(values = c("Yes" = 0.9, "No" = 0.2)) +
  geom_label_repel(data = mr_volcano_dat[(p<0.05 & as.character(bmi_metabolite)=="Yes") | p < 5e-7], aes(label = label), size = 2.5,
                   color="black", max.overlaps = 6, show.legend = F) +
  theme_classic() +
  theme(legend.position = "top") +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    alpha = guide_legend(override.aes = list(size = 3))
  ) +
  labs(x = "MR estimate - OR for all-cause HF per SD metabolite",
       y = expression(-log[10](p.value)),
       color = "BMI metabolite",
       alpha = "BMI metabolite")
p

pathway_colors <- c("Lipid metabolism" = "#E41A1C",
                    "Lipoproteins" ="#377EB8",
                    "Amino acid metabolism" = "#4DAF4A",
                    "Carbohydrate/energy metabolism" = "#984EA3",
                    "Nucleotide metabolism" = "#999999",
                    "Cofactors and vitamins" = "#FF7F00",
                    "Inflammation" = "#FFFF33",
                    "Fluid balance" = "#F781BF",
                    "Xenobiotics/other" = "#A65628",
                    "Unknown")

lims_x <- c(0.5, 2)
p_list <- lapply(split(forest_dat[!grepl("^X", as.character(exposure_lab)) & p < 0.05], by="outcome"),
                 function(x) {
  ggplot(x, aes(y = exposure_lab, x = or, color=pathway_group)) +
    geom_vline(xintercept = 1, color="darkgray") +
    geom_errorbarh(aes(xmin = ifelse(or_lb<lims_x[1], lims_x[1], or_lb),
                       xmax = ifelse(or_ub>lims_x[2], lims_x[2], or_ub)), color="lightgray", height=0) +
    geom_point(aes(shape = p_sig), size=2) +
    scale_color_manual(values = pathway_colors) +
    scale_shape_manual(values = c("P < 0.05"=19, "FDR P < 0.05" = 17)) +
    labs(x = "OR for outcome per 1 SD increase in metabolite",
         shape = "P-value",
         size = "P-value",
         color = "Pathway") +
    lims(x = lims_x) +
    theme_light() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size=6),
      axis.title.y = element_blank()
    ) +
    facet_wrap(~outcome, ncol=1, scales = "free_y")
})

plot_heights <- sapply(p_list, function(p) {
  built <- ggplot_build(p)
  y_counts <- sapply(built$data, function(layer) {
    if (!is.null(layer$y)) length(unique(layer$y)) else 0
  })
  max(y_counts)
})
plot_heights <- plot_heights / sum(plot_heights)
base_strip_height <- 0.20  # tweak this
plot_heights_adj <- plot_heights * 1 + base_strip_height
plot_heights_adj <- plot_heights_adj / sum(plot_heights_adj)



combined_legend <- get_legend(
  ggplot(forest_dat[!grepl("^X", as.character(exposure_lab)) & p < 0.05],
         aes(y = exposure_lab, x = or, color=pathway_group)) +
    geom_point(aes(shape = p_sig)) +
    scale_color_manual(values = pathway_colors) + #palette = "Set1", direction = 1) +
    scale_shape_manual(values = c("P < 0.05"=19, "FDR P < 0.05" = 17)) +
    labs(x = "OR for outcome per 1 SD increase in metabolite",
         shape = "P-value",
         size = "P-value",
         color = "Pathway") +
    theme_light() +
    guides(
      color = guide_legend(override.aes = list(size = 3)),
      shape = guide_legend(override.aes = list(size = 3))
    ) +
    facet_wrap(~outcome)
)


forest_combined <- plot_grid(plotlist = lapply(seq_along(p_list), function(i) {
                                p <- p_list[[i]] + theme(legend.position = "none")
                                if (i < 3) p <- p + theme(axis.title.x = element_blank())
                                p
                             }),
                             ncol = 1,
                             rel_heights = plot_heights_adj,
                             align = "v",
                             axis = "lr")
forest_legend <- plot_grid(forest_combined, combined_legend, nrow=1, rel_widths = c(1,0.3))
p0 <- plot_grid(p, forest_legend, nrow = 2, rel_heights = c(0.8, 1))
p0

# save
ggsave(outcome_volcano_forest, p0, width = 10, height = 14, dpi=300, bg="white")



# IVW MR significant metabolites
significant_dat <- copy(outcomes)
significant_dat[, ivw_pos := p[method=="mr_ivw"] < 0.05, by=c("label", "platform", "outcome")]
significant_dat[is.na(b), `:=`(b=0, b_se=0, missing=TRUE)]
significant_dat[is.na(missing), `:=`(missing=FALSE)]
significant_dat[, sig := p < 0.05]
significant_dat[, fill_col := ifelse(sig, as.character(pathway_group), NA)]
significant_dat <- significant_dat[, .SD[any(ivw_pos, na.rm=TRUE)], by=c("label", "platform", "outcome")]
significant_dat[, `:=`(or   = exp(b),
                       orlb = exp(b-1.96*b_se),
                       orub = exp(b+1.96*b_se))]

significant_dat[
  , exposure_lab := factor(
    exposure_lab,
    levels = rev(unique(exposure_lab[order(pathway_group, -b)]))
  )
]

# report IVW sig and MR consistent
significant_dat[, .(num=uniqueN(exposure_lab))]
significant_dat[, .(num=uniqueN(exposure_lab)), by=c("outcome")][order(outcome)]
significant_dat[, .(num=uniqueN(exposure_lab)), by=c("outcome","pathway_group")][order(outcome, pathway_group)]



lower_lim <- 0.5
upper_lim <- 1.5
dodge <- 0.8

p <- ggplot(significant_dat,
            aes(y = exposure_lab,
                x = or,
                color=pathway_group,
                shape=method,
                fill = factor(fill_col, levels = levels(pathway_group)),
                alpha=missing)) +
  geom_vline(xintercept = 1, color="darkgray") +
  geom_errorbarh(aes(xmin = ifelse(orlb<lower_lim,lower_lim,orlb), xmax = ifelse(orub>upper_lim,upper_lim,orub)),
                 height=0,
                 color="lightgray",
                 position = position_dodge(width=dodge)) +
  geom_point(size=2, position = position_dodge(width=dodge)) +
  scale_shape_manual(values = c(
    "mr_ivw" = 21,
    "mr_egger" = 22,
    "mr_weighted_median" = 24,
    "mr_weighted_mode" = 25
  )) +
  scale_fill_brewer(palette = "Set1", na.value = "white") +
  scale_color_brewer(palette = "Set1") +
  scale_alpha_manual(values = c("TRUE"=0, "FALSE"=1)) +
  scale_y_discrete(labels = function(x) {
    clean <- gsub("_", " ", sub("^mr_","", x))
    ifelse(clean=="ivw", "IVW", stringr::str_to_title(clean))
  }) +
  guides(fill = "none", alpha="none") +
  labs(x = "OR per 1 SD increase in Metabolite for HFpEF", color = "Pathway", fill = "Pathway") +
  theme_minimal(base_size = 8) +
  lims(x=c(lower_lim,upper_lim)) +
  facet_wrap(~outcome, scales = "free_x") +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 5),
    axis.title.y = element_blank()
  )
p

ggsave(outcome_mr_sig, p, width = 9, height = 10, dpi=300, bg="white")



# look at loci overlap
metab_loci <- fread(metab_loci_fp, select = c("trait", "id", "rsid", "chr", "bp", "ea", "oa", "eaf", "beta", "se", "p"))
metab_loci[map, `:=`(label = i.label, pathway_group = i.pathway_group, platform = i.platform), on=c("trait"="gwas_name", "id"="gwas_id")]
metab_loci[, pathway_group := factor(pathway_group, levels = c("Lipid metabolism",
                                                                 "Lipoproteins",
                                                                 "Amino acid metabolism",
                                                                 "Carbohydrate/energy metabolism",
                                                                 "Nucleotide metabolism",
                                                                 "Cofactors and vitamins",
                                                                 "Inflammation",
                                                                 "Fluid balance",
                                                                 "Xenobiotics/other",
                                                                 "Unknown"))]
metab_loci[rsid=="", rsid := NA_character_]
metab_loci <- metab_loci[!is.na(label)] # only compare metabolites that the studies used - i.e. those in the mapping
metab_loci <- unique(metab_loci)
metab_loci[, trait := paste0(label, " [", platform, "]")]
metab_loci[, bmi_metabolite := ifelse(trait %in% rep[replicates_trials_mr=="Trial & MR consistent", label_plat], TRUE, FALSE)]

# save the BMI-metabolite variants
metab_loci[, `:=`(chr_bp_ea_oa = paste0(chr, ":", bp, "_", pmin(ea, oa), "_", pmax(ea, oa)),
                  new_ea = pmin(ea, oa))]
metab_loci[ea!=new_ea, `:=`(beta = -1 * beta, eaf = 1 - eaf)]
bmi_metabolite_variants <- metab_loci[rsid!="" & !is.na(rsid),
                                      .(num_traits = uniqueN(trait),
                                        metabolite_traits = paste0(paste0(trait, "#gwas_id=", id, "#eaf=", eaf, "#beta=", beta, "#se=", se, "#p=", p), collapse=";;"),
                                        category = fcase(all(bmi_metabolite), "bmi_metabolite_variant",
                                                         all(!bmi_metabolite), "non_bmi_metabolite_variant",
                                                         default = "both")), by = c("rsid", "chr_bp_ea_oa")]

fwrite(bmi_metabolite_variants[order(category)], bmi_metabolite_snps, sep="\t")


library(ggvenn)
venn_list <- list(
  `Non BMI Metabolite` = unique(metab_loci[bmi_metabolite == FALSE, rsid]),
  `BMI Metabolite` = unique(metab_loci[bmi_metabolite == TRUE, rsid])
)
p <- ggvenn(venn_list,
       fill_color = c("#E69F00", "#56B4E9"),
       stroke_size = 0.5,
       set_name_size = 2.8,
       auto_scale = TRUE)
p

ggsave(bmi_metabolite_variant_venn, p, width = 6, height = 5, dpi=300, bg="white")


# traits to process (all BMI variants with an instrument)
bmi_metabolites_with_mr <- outcomes[!is.na(b), unique(paste0(label, " [", platform, "]"))]
sig_trait_snps <- lapply(bmi_metabolites_with_mr, function(t) {

  rsids      <- metab_loci[trait == t, rsid]
  trait_snps <- metab_loci[rsid %in% rsids]
  trait_snps[, .(sig_trait = t, sig_trait_nsnp = length(rsids), overlap_trait = trait, overlap_trait_pathway = pathway_group, rsid)]

}) |> rbindlist()



overlap_counts <- sig_trait_snps[, .(
  n_shared_snps = uniqueN(rsid)
), by = .(sig_trait, overlap_trait)]

sig_trait_nsnp <- sig_trait_snps[, .(sig_trait_nsnp = uniqueN(rsid)), by = .(sig_trait)]
pct_overlap <- merge(overlap_counts, sig_trait_nsnp, by = c("sig_trait"))
pct_overlap[, pct_overlap_snps := n_shared_snps / sig_trait_nsnp]

pct_overlap[, sig_trait := factor(sig_trait, levels = unique(sig_trait))]
pct_overlap[, overlap_trait := factor(overlap_trait, levels = unique(metab_loci$trait))]
pct_overlap <- tidyr::complete(pct_overlap, sig_trait, overlap_trait, fill = list(pct_overlap_snps=0)) |> as.data.table()
pct_overlap[metab_loci, overlap_pathway_group := i.pathway_group, on=c("overlap_trait"="trait")]
pct_overlap[metab_loci, sig_trait_pathway_group := i.pathway_group, on=c("sig_trait"="trait")]
pct_overlap[sig_trait_snps, nsnp := i.sig_trait_nsnp, on=c("sig_trait")]
pct_overlap[, overlap_trait := factor(overlap_trait, levels = unique(overlap_trait[order(overlap_pathway_group)]))]
pct_overlap[, sig_trait_lab := paste0(sig_trait, "[nsnp=", nsnp, "]")]
pct_overlap[, overlap_sum := mean(pct_overlap_snps, na.rm=T), by="sig_trait_lab"]
pct_overlap[, sig_trait_lab := factor(sig_trait_lab, levels = unique(sig_trait_lab[order(-overlap_sum)]))]


library(RColorBrewer)
dir.create(outcome_snp_overlap_dir, showWarnings = FALSE, recursive = TRUE)
for (ii in as.character(unique(pct_overlap$sig_trait_lab))) {
  dat <- copy(pct_overlap)
  dat <- dat[as.character(sig_trait_lab)==ii & as.character(sig_trait) != as.character(overlap_trait), ]
  p <- ggplot(dat,
              aes(x = overlap_trait,
                  y = pct_overlap_snps,
                  group = overlap_pathway_group,
                  color=overlap_pathway_group,
                  fill = overlap_pathway_group)
  ) +
    geom_area(position = "identity", alpha = 0.6, color = NA) +
    geom_line(aes(group = sig_trait), linewidth = 0.3) +
    scale_fill_discrete(drop = FALSE) +
    scale_color_discrete(drop = FALSE)+
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(1, "lines"),
      panel.grid = element_blank(),
      strip.text = element_text(size=5)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      subtitle = ii,
      x = paste0("All Metabolites with instruments (n=", metab_loci[, uniqueN(id)], ")"),
      y = "Proportion of Overlapping SNPs",
      color = "Pathway",
      fill = "Pathway"
    )

  ii_clean <- gsub("[/\\:*?\"<>|]", "_", ii)
  fp <- file.path(outcome_snp_overlap_dir, paste0(ii_clean, "_snp_overlap.png"))
  ggsave(fp, p, width = 10, height = 5, dpi=300, bg="white")
}


outcome_instr_summ <- pct_overlap[!is.na(sig_trait),
  {
    oth_traits <- overlap_trait[as.character(overlap_trait) != as.character(sig_trait) & !is.na(n_shared_snps)]
    oth_pct    <- pct_overlap_snps[as.character(overlap_trait) != as.character(sig_trait) & !is.na(n_shared_snps)]

    if (length(oth_traits)) {
      o <- order(-oth_pct)
      overlap_str <- paste0(oth_traits[o], " (", round(oth_pct[o]*100), "%)", collapse = "; ")
    } else {
      overlap_str <- NA_character_
    }

    .(
      instrument_size    = nsnp[1],
      num_overlap_traits = length(oth_traits),
      overlap_traits     = overlap_str,
      max_overlap_snps_pct   = ifelse(is.infinite(max(oth_pct)), NA_real_, max(oth_pct)),
      mean_overlap_snps_pct  = ifelse(is.nan(mean(oth_pct)), NA_real_, mean(oth_pct))
    )
  },
  by = .(sig_trait)
]
outcome_instr_summ[metab_loci, `:=`(pathway_group = i.pathway_group), on=c("sig_trait"="trait")]


# save
fwrite(outcome_instr_summ[order(num_overlap_traits),
                            .(metabolite_trait      = sig_trait,
                              metabolite_trait_nsnp = instrument_size,
                              pathway_group         = pathway_group,
                              n_overlap_metabolites = num_overlap_traits,
                              overlap_metabolites   = overlap_traits
                              )], outcome_instr_tbl, sep="\t")

# pct overlap by pathway
pathway_overlap_formatted <- pct_overlap[as.character(sig_trait) != as.character(overlap_trait), .(
  overlap_pathway_group = "all other",
  min_overlap_pct = sprintf("%.1f%%", 100*min(pct_overlap_snps, na.rm=T)),
  max_overlap_pct = sprintf("%.1f%%", 100*max(pct_overlap_snps, na.rm=T)),
  median_overlap_pct = sprintf("%.1f%%", 100*median(pct_overlap_snps, na.rm=T)),
  mean_overlap_pct = sprintf("%.1f%%", 100*mean(pct_overlap_snps, na.rm=T))
), by = c("sig_trait_pathway_group")]

# Format pathway_pathway_overlap
pathway_pathway_overlap_formatted <- pct_overlap[as.character(sig_trait) != as.character(overlap_trait), .(
  min_overlap_pct = sprintf("%.1f%%", 100*min(pct_overlap_snps, na.rm=T)),
  max_overlap_pct = sprintf("%.1f%%", 100*max(pct_overlap_snps, na.rm=T)),
  median_overlap_pct = sprintf("%.1f%%", 100*median(pct_overlap_snps, na.rm=T)),
  mean_overlap_pct = sprintf("%.1f%%", 100*mean(pct_overlap_snps, na.rm=T))
), by = c("sig_trait_pathway_group", "overlap_pathway_group")]

fwrite(rbind(pathway_overlap_formatted, pathway_pathway_overlap_formatted), pathway_overlap_tbl, sep="\t")





# Overlap with Lui et al.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10942215/

hf_types <- c("all", "hfref", "hfpef")
models <- 1:3
cols <- c("hr", "ci_95", "p_fdr")
vector <- as.vector(unlist(lapply(hf_types, function(hf) {sapply(models, function(m) {paste0(hf, "_model_", m, "_", cols)})})))

# MS
liu_fp <- "/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/.results/tables/outcomes/circhf_circhf-2023-010896_supp2.xlsx"
hf1 <- read_xlsx(liu_fp, sheet=1, skip=4, col_names = c("liu_name", "pathway", "subpathway", "hmdb", vector, "new", "fig_num", "lasso")) |> as.data.table()
hf2 <- read_xlsx(liu_fp, sheet=2, skip=4, col_names = c("liu_name", vector))|> as.data.table()
hf_obs  <- rbind(hf1, hf2, fill=TRUE)
hf_obs[map, `:=`(label = i.label, platform = "metabolon"), on="liu_name"]
hf_obs[is.na(label), `:=`(label = liu_name, platform = "metabolon")]
setcolorder(hf_obs, c("label", "platform"))


# NMR
julkunen_fp <- "/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/.results/tables/outcomes/julkunen_2023_summary_statistics.csv"
hf_nmr <- fread(julkunen_fp)
hf_nmr[map, `:=`(label = i.label, platform = "nightingale"), on=c("biomarker_name"="julkunen_name")]

# join
hf_obs <- rbind(
  hf_obs,
  hf_nmr[, .(label, platform, julkunen_name = biomarker_name, pathway = group_name,
             all_model_1_hr = hazard_ratio,
             all_model_1_ci_95 = ifelse(is.na(ci_lower), NA_character_, sprintf("%.2f, %.2f", ci_lower, ci_upper)),
             all_model_1_p_fdr = p.adjust(pvalue, method="fdr"), # reported P is unadjusted
             hfref_model_1_hr = 1,
             hfref_model_1_ci_95 = "",
             hfref_model_1_p_fdr = 1,
             hfpef_model_1_hr = 1,
             hfpef_model_1_ci_95 = "",
             hfpef_model_1_p_fdr = 1)],
  fill = TRUE
)


bmi_metabs <- merge(rep[replicates_trials_mr=="Trial & MR consistent",
                        .(label, platform, label_plat,
                          bmi_estimate_bbs = as.numeric(as.numeric(bmi_estimate_bbs)),
                          bmi_estimate_direct = as.numeric(as.numeric(bmi_estimate_direct)),
                          bmi_estimate_mr = as.numeric(as.numeric(bmi_estimate_mr)),
                          bmi_std.error_bbs = as.numeric(as.numeric(bmi_std.error_bbs)),
                          bmi_std.error_direct = as.numeric(as.numeric(bmi_std.error_direct)),
                          bmi_std.error_mr = as.numeric(as.numeric(bmi_std.error_mr)),
                          bmi_p.value_fdr_bbs = as.numeric(as.numeric(bmi_p.value_fdr_bbs)),
                          bmi_p.value_interaction_fdr_direct = as.numeric(as.numeric(bmi_p.value_interaction_fdr_direct)),
                          bmi_p.value_fdr_mr = as.numeric(as.numeric(bmi_p.value_fdr_mr))
                          )],
                    hf_obs[, .SD, .SDcols = c("label", "platform", grep("hr|p_fdr|ci_95", names(hf_obs), value=T))],
                    by = c("label", "platform"),
                    all.x = TRUE)
bmi_metabs <- merge(bmi_metabs,
                    outcomes[method=="mr_ivw" & outcome=="All-cause HF", .(label, platform, estimate_hfall_mr=b, std.error_hfall_mr=b_se, p.value_hfall_mr=p)],
                    by = c("label", "platform"),
                    all.x = TRUE)
bmi_metabs <- merge(bmi_metabs,
                    outcomes[method=="mr_ivw" & outcome=="Non-ischaemic HFrEF", .(label, platform, estimate_hfref_mr=b, std.error_hfref_mr=b_se, p.value_hfref_mr=p)],
                    by = c("label", "platform"),
                    all.x = TRUE)
bmi_metabs <- merge(bmi_metabs,
                    outcomes[method=="mr_ivw" & outcome=="Non-ischaemic HFpEF", .(label, platform, estimate_hfpef_mr=b, std.error_hfpef_mr=b_se, p.value_hfpef_mr=p)],
                    by = c("label", "platform"),
                    all.x = TRUE)

fwrite(bmi_metabs[, .(label, platform,
                      bmi_estimate_bbs, bmi_std.error_bbs, bmi_p.value_fdr_bbs,
                      bmi_estimate_direct, bmi_std.error_direct, bmi_p.value_interaction_fdr_direct,
                      bmi_estimate_mr, bmi_std.error_mr, bmi_p.value_fdr_mr,
                      observational_source = fcase(!is.na(all_model_1_hr) & platform=="nightingale", "Julkunen et al. 2023",
                                                   !is.na(all_model_1_hr) & platform=="metabolon",   "Liu et al. 2024",
                                                   default = NA_character_),
                      obs_hf_all_hr = all_model_1_hr,
                      obs_hf_all_95ci = all_model_1_ci_95,
                      obs_hf_all_p.value_fdr = all_model_1_p_fdr,
                      mr_hf_all_or      = exp(estimate_hfall_mr),
                      mr_hf_all_95ci    = ifelse(is.na(estimate_hfall_mr), NA_character_, sprintf("%.2f, %.2f", exp(estimate_hfall_mr - 1.96*std.error_hfall_mr), exp(estimate_hfall_mr + 1.96*std.error_hfall_mr))),
                      mr_hf_all_p.value = p.value_hfall_mr,
                      obs_hfref_hr = hfref_model_1_hr,
                      obs_hfref_95ci = hfref_model_1_ci_95,
                      obs_hfref_p.value_fdr = hfref_model_1_p_fdr,
                      mr_hfref_or       = exp(estimate_hfref_mr),
                      mr_hfref_95ci    = ifelse(is.na(estimate_hfall_mr), NA_character_, sprintf("%.2f, %.2f", exp(estimate_hfref_mr - 1.96*std.error_hfref_mr), exp(estimate_hfref_mr + 1.96*std.error_hfref_mr))),
                      mr_hfref_p.value = p.value_hfref_mr,
                      obs_hfpef_hr = hfpef_model_1_hr,
                      obs_hfpef_95ci = hfpef_model_1_ci_95,
                      obs_hfpef_p.value_fdr = hfpef_model_1_p_fdr,
                      mr_hfpef_or       = exp(estimate_hfpef_mr),
                      mr_hfpef_95ci    = ifelse(is.na(estimate_hfall_mr), NA_character_, sprintf("%.2f, %.2f", exp(estimate_hfpef_mr - 1.96*std.error_hfpef_mr), exp(estimate_hfpef_mr + 1.96*std.error_hfpef_mr))),
                      mr_hfpef_p.value = p.value_hfpef_mr)],
       "/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/.results/tables/outcomes/observational_hf_overlap.tsv",
       sep="\t")

# get overlap
library(eulerr)
bmi_metabs_euler_dat <- merge(rep[, .(label, platform, replicates_trials_mr)],
                              hf_obs[, .SD, .SDcols = c("label", "platform", grep("hr|p_fdr|ci_95", names(hf_obs), value=T))],
                              by  = c("label", "platform"),
                              all = TRUE)
bmi_metabs_euler_dat <- bmi_metabs_euler_dat[, {
  .(
    label            = label,
    All              = TRUE,
    `MS platform`    = platform=="metabolon",
    `NMR platform`   = platform=="nightingale",
    `BMI-metabolite` = ifelse(is.na(replicates_trials_mr),FALSE, replicates_trials_mr=="Trial & MR consistent"),
    `All-cause HF`   = ifelse(is.na(all_model_1_p_fdr),   FALSE, all_model_1_p_fdr   < 0.05),
    HFpEF            = ifelse(is.na(hfpef_model_1_p_fdr), FALSE, hfpef_model_1_p_fdr < 0.05),
    HFrEF            = ifelse(is.na(hfref_model_1_p_fdr), FALSE, hfref_model_1_p_fdr < 0.05)
  )}]

#print out
bmi_metabs_euler_dat[, c(lapply(.SD, sum), lapply(.SD, function(x) sum(x)/.N)), by="BMI-metabolite", .SDcols = setdiff(names(bmi_metabs_euler_dat), "label")]

set.seed(122)
fit <- euler(bmi_metabs_euler_dat[`MS platform`==TRUE, .SD, .SDcols = setdiff(names(bmi_metabs_euler_dat), c("label", "All", "NMR platform"))],
             input = "disjoint",
             shape = "ellipse")
p <- plot(fit,
          fills = list(fill = c("white","#E69F00", "darkgreen", "royalblue", "darkred"), alpha = 0.3),
          legend = list(alpha=1),
          labels = list(fontsize = 7),
          quantities = TRUE,
          edges = TRUE)
p

png(liu_model1_overlap_venn_ms,
    width = 6, height = 5, units = "in", res = 300)
grid::grid.draw(p)
dev.off()

set.seed(125)
fit <- euler(bmi_metabs_euler_dat[`NMR platform`==TRUE, .SD, .SDcols = setdiff(names(bmi_metabs_euler_dat), c("label", "All", "MS platform"))],
             input = "disjoint",
             shape = "ellipse")
p <- plot(fit,
          fills = list(fill = c("white","#E69F00", "darkgreen", "royalblue", "darkred"), alpha = 0.3),
          legend = list(alpha=1),
          labels = list(fontsize = 7),
          quantities = TRUE,
          edges = TRUE)
p

png(liu_model1_overlap_venn_nmr,
    width = 6, height = 5, units = "in", res = 300)
grid::grid.draw(p)
dev.off()



# plot concordance
sig_colors <- list(
  pos_sig  = "royalblue",
  pos_nsig = adjustcolor("royalblue", alpha.f = 0.5),
  neg_sig  = "firebrick",
  neg_nsig = adjustcolor("firebrick", alpha.f = 0.5),
  default  = adjustcolor("gray", alpha.f = 0.5)
)

bmi_metabs_long <- bmi_metabs[, .(
  label_plat,
  bbs = fcase(bmi_estimate_bbs > 0 & bmi_p.value_fdr_bbs <  0.05, sig_colors[["pos_sig"]],
              bmi_estimate_bbs > 0 & bmi_p.value_fdr_bbs >= 0.05, sig_colors[["pos_nsig"]],
              bmi_estimate_bbs < 0 & bmi_p.value_fdr_bbs <  0.05, sig_colors[["neg_sig"]],
              bmi_estimate_bbs < 0 & bmi_p.value_fdr_bbs >= 0.05, sig_colors[["neg_nsig"]],
              default = sig_colors[["default"]]),
  dir = fcase(bmi_estimate_direct > 0 & bmi_p.value_interaction_fdr_direct <  0.05, sig_colors[["pos_sig"]],
              bmi_estimate_direct > 0 & bmi_p.value_interaction_fdr_direct >= 0.05, sig_colors[["pos_nsig"]],
              bmi_estimate_direct < 0 & bmi_p.value_interaction_fdr_direct <  0.05, sig_colors[["neg_sig"]],
              bmi_estimate_direct < 0 & bmi_p.value_interaction_fdr_direct >= 0.05, sig_colors[["neg_nsig"]],
              default = sig_colors[["default"]]),
  bmi = fcase(bmi_estimate_mr > 0 & bmi_p.value_fdr_mr <  0.05, sig_colors[["pos_sig"]],
              bmi_estimate_mr > 0 & bmi_p.value_fdr_mr >= 0.05, sig_colors[["pos_nsig"]],
              bmi_estimate_mr < 0 & bmi_p.value_fdr_mr <  0.05, sig_colors[["neg_sig"]],
              bmi_estimate_mr < 0 & bmi_p.value_fdr_mr >= 0.05, sig_colors[["neg_nsig"]],
              default = sig_colors[["default"]]),
  mr_all = fcase(estimate_hfall_mr > 0 & p.value_hfall_mr <  0.05, sig_colors[["pos_sig"]],
                 estimate_hfall_mr > 0 & p.value_hfall_mr >= 0.05, sig_colors[["pos_nsig"]],
                 estimate_hfall_mr < 0 & p.value_hfall_mr <  0.05, sig_colors[["neg_sig"]],
                 estimate_hfall_mr < 0 & p.value_hfall_mr >= 0.05, sig_colors[["neg_nsig"]],
                 default = sig_colors[["default"]]),
  mr_hfref = fcase(estimate_hfref_mr > 0 & p.value_hfref_mr <  0.05, sig_colors[["pos_sig"]],
                   estimate_hfref_mr > 0 & p.value_hfref_mr >= 0.05, sig_colors[["pos_nsig"]],
                   estimate_hfref_mr < 0 & p.value_hfref_mr <  0.05, sig_colors[["neg_sig"]],
                   estimate_hfref_mr < 0 & p.value_hfref_mr >= 0.05, sig_colors[["neg_nsig"]],
                   default = sig_colors[["default"]]),
  mr_hfpef = fcase(estimate_hfpef_mr > 0 & p.value_hfpef_mr <  0.05, sig_colors[["pos_sig"]],
                   estimate_hfpef_mr > 0 & p.value_hfpef_mr >= 0.05, sig_colors[["pos_nsig"]],
                   estimate_hfpef_mr < 0 & p.value_hfpef_mr <  0.05, sig_colors[["neg_sig"]],
                   estimate_hfpef_mr < 0 & p.value_hfpef_mr >= 0.05, sig_colors[["neg_nsig"]],
                   default = sig_colors[["default"]]),
  obs_all = fcase(all_model_1_hr > 1 & all_model_1_p_fdr <  0.05, sig_colors[["pos_sig"]],
                  all_model_1_hr > 1 & all_model_1_p_fdr >= 0.05, sig_colors[["pos_nsig"]],
                  all_model_1_hr < 1 & all_model_1_p_fdr <  0.05, sig_colors[["neg_sig"]],
                  all_model_1_hr < 1 & all_model_1_p_fdr >= 0.05, sig_colors[["neg_nsig"]],
                  default = sig_colors[["default"]]),
  obs_hfref = fcase(hfref_model_1_hr > 1 & hfref_model_1_p_fdr <  0.05, sig_colors[["pos_sig"]],
                    hfref_model_1_hr > 1 & hfref_model_1_p_fdr >= 0.05, sig_colors[["pos_nsig"]],
                    hfref_model_1_hr < 1 & hfref_model_1_p_fdr <  0.05, sig_colors[["neg_sig"]],
                    hfref_model_1_hr < 1 & hfref_model_1_p_fdr >= 0.05, sig_colors[["neg_nsig"]],
                    default = sig_colors[["default"]]),
  obs_hfpef = fcase(hfpef_model_1_hr > 1 & hfpef_model_1_p_fdr <  0.05, sig_colors[["pos_sig"]],
                    hfpef_model_1_hr > 1 & hfpef_model_1_p_fdr >= 0.05, sig_colors[["pos_nsig"]],
                    hfpef_model_1_hr < 1 & hfpef_model_1_p_fdr <  0.05, sig_colors[["neg_sig"]],
                    hfpef_model_1_hr < 1 & hfpef_model_1_p_fdr >= 0.05, sig_colors[["neg_nsig"]],
                    default = sig_colors[["default"]])
)]
bmi_metabs_long <- melt(bmi_metabs_long, id.vars = "label_plat")
bmi_metabs_long[, facet := fcase(grepl("bbs|dir|bmi", variable), "base",
                                 grepl("all", variable), "all",
                                 grepl("hfref", variable), "hfref",
                                 grepl("hfpef", variable), "hfpef")]

bmi_metabs_long[, variable := factor(variable, levels = c("bbs","dir","bmi","obs_all","mr_all","obs_hfpef","mr_hfpef","obs_hfref","mr_hfref"))]
bmi_metabs_long[, sum := sum(match(value, sig_colors)), by = c("label_plat", "facet")]
bmi_metabs_long[, label_plat := factor(label_plat, levels = unique(label_plat[facet=="base"][order(sum[facet=="base"])]))]
bmi_metabs_long[, facet := factor(facet, levels = c("base", "all", "hfpef", "hfref"), labels = c("BMI-metabolite\nselection", "Outcome\nAll-cause HF", "Outcome\nHFpEF", "Outcome\nHFrEF"))]

p <- ggplot(bmi_metabs_long, #[!grepl("^X", label_plat)],
       aes(x = variable, y = label_plat, fill = value)) +
  geom_tile(color="white") +
  scale_fill_identity(
    name = "Direction & Sig.",
    breaks = sig_colors,
    labels = c("\u03B2+ & P<0.05", "\u03B2+ & P\u22650.05", "\u03B2- & P<0.05", "\u03B2- & P\u22650.05", "Missing"),
    guide = "legend"
  ) +
  scale_x_discrete(labels = c(bbs = "BMI in bariatric surgery",
                              dir = "BMI in calorie restriction",
                              bmi = "Lifetime BMI MR",
                              mr_all    = "MR",
                              mr_hfref  = "MR",
                              mr_hfpef  = "MR",
                              obs_all   = "Observational",
                              obs_hfref = "Observational",
                              obs_hfpef = "Observational")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title  = element_blank()
  ) +
  facet_wrap(~facet, nrow=1, scales="free_x")
p

ggsave(out_mr_obs_heatmap,
       p, width = 8, height = 14, dpi=300, bg="white")














#
#
#
# # R2 of sig MR results
# sig_mr_exposure_snps <- trait_loci[unique(hf[p<0.05, .(exposure)]), on = c("trait"="exposure"), nomatch=NULL]
#
# # write out the variants file
# variants_file <- tempfile()
# data.table::fwrite(data.table::data.table(RSID=unique(metab_loci$rsid)), variants_file, sep="\t", col.names=FALSE)
#
# sig_instr_snps_file <- tempfile()
# data.table::fwrite(data.table::data.table(RSID=unique(sig_mr_exposure_snps$rsid)), sig_instr_snps_file, sep="\t", col.names=FALSE)
#
#
# # create and run the LD matrix plink2 function
# plink_ldmat <- tempfile()
# cmd <- paste("/usr/local/plink/plink2",
#              "--pfile", "/Users/xx20081/Documents/local_data/genome_reference/ukb_reference_genome/uk10k",
#              "--extract", variants_file,
#              "--ld-snp-list", sig_instr_snps_file,
#              "--keep-allele-order ",
#              "--r2-phased", "inter-chr",
#              #"--ld-window-kb", 9999999,
#              "--ld-window-r2", -1,
#              "--out", plink_ldmat)
# system(cmd)
#
# # read in the LD matrix variant rsids
# ld_dat <- data.table::fread(paste0(plink_ldmat, ".vcor"), sep="\t")
#
# # clean up
# unlink(paste0(plink_ldmat, ".vcor"))
# unlink(variants_file)
# unlink(sig_instr_snps_file)
#
# # join to traits
# ld_dat <- ld_dat[, .(sig_instr_snp = ID_A, other_instr_snp = ID_B, r2 = PHASED_R2)]
#
#
# AB <- merge(sig_mr_exposure_snps[, .(key=1, rsid, trait, platform)],
#             metab_loci[, .(key=1, rsid, trait, platform)],
#             by = "key", allow.cartesian = TRUE)[, key := NULL]
#
# AB[ld_dat, r2 := r2, on=c("rsid.x"="sig_instr_snp", "rsid.y"="other_instr_snp")]
# AB[map, pathway.x := i.pathway, on=c("trait.x"="gwas_name", "platform.x"="platform")]
# AB[map, pathway.y := i.pathway, on=c("trait.y"="gwas_name", "platform.y"="platform")]
#
# avg_r2 <- AB[, .(n_snps.x = uniqueN(rsid.x),
#                  n_snps.y = uniqueN(rsid.y),
#                  high_r2  = sum(r2 > 0.5, na.rm=T)#, # number of SNP pairs with r2>0.1
#                  #max_r2   = max(r2, na.rm=T),
#                  #mean_r2  = mean(r2, na.rm=T)
# ), by=c("trait.x","platform.x","trait.y","platform.y", "pathway.y")]
# avg_r2[, trait.x_lab := paste(trait.x, platform.x, sep="|")]
# avg_r2[, trait.y_lab := paste(trait.y, platform.y, sep="|")]
#
# avg_r2[, num_sig_instr_highr2 := sum(high_r2, na.rm=T), by="trait.x_lab"]
# avg_r2[, num_other_instr_highr2 := sum(high_r2, na.rm=T), by="trait.y_lab"]
# avg_r2[, trait.x_lab := factor(trait.x_lab, levels = unique(trait.x_lab[order(-num_sig_instr_highr2)]))]
# avg_r2[, trait.y_lab := factor(trait.y_lab, levels = unique(trait.y_lab[order(pathway.y)]))]
#
# ggplot(avg_r2, aes(y = trait.x_lab, x = trait.y_lab, fill = pathway.y, color = high_r2>0)) +
#   geom_tile(alpha=0.8) +
#   # geom_point(data = avg_r2[high_r2>0], color="red", size=1) +
#   scale_fill_viridis_d() +
#   scale_color_manual(values = c("TRUE"="red", "FALSE"=NA)) +
#   theme(
#     axis.text.x = element_blank() # element_text(angle=45, hjust=1)
#   ) +
#   labs(
#     x = "Metabolites",
#     y = "Weight loss trial & BMI MR & Heart failure outcome\nsignificant metabolites",
#     color = "Any metabolite instrument SNP pair R2 > 0.1",
#     fill = "Pathway"
#   )
#
#
# bar_dat <- avg_r2[, .(other_instr_has_highr2 = sum(high_r2 > 0, na.rm=T)), by=c("trait.x_lab", "trait.y_lab", "pathway.y")]
# bar_dat <- bar_dat[, .(num_other_instr_has_highr2 = sum(other_instr_has_highr2)), by=c("trait.x_lab", "pathway.y")]
#
# ggplot(bar_dat, aes(y = trait.x_lab, x = num_other_instr_has_highr2, fill = pathway.y)) +
#   geom_col() +
#   labs(y = "Weight loss & BMI MR replicating & HF relevant metabolites",
#        x = "Number of total metabolites with at least 1 SNP in LD with instrument (R2>0.1)",
#        fill = "Pathway")
#
#
#
#
#
#
#
#
#
# # look at overlap by significant MR
# for (out in unique(hf$outcome)) {
#   trait_loci[, (out) := ifelse(trait %in% hf[outcome==out & p<0.05, exposure], TRUE, FALSE)]
# }
# by_outcome_trait_loci <- melt(trait_loci, id.vars = c("trait", "overall_clump_lab", "num_traits"), measure.vars = as.character(unique(hf$outcome)))[value==TRUE]
#
# ggplot(by_outcome_trait_loci, aes(x = overall_clump_lab, y = trait, fill = num_traits)) +
#   geom_tile() +
#   scale_fill_viridis() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust=1)
#   ) +
#   facet_wrap(~variable, scales = "free_y")
#
#
#
#
#
#
# # compare MR results between phenotypes
# mr_comp <- dcast(hf, exposure + metabolite_platform ~ outcome, value.var = c("b", "b_se"))
#
# xout        <- "af"
# yout        <- "dcm_mtag"
# error_ratio <- mean(mr_comp[, get(paste0("b_se_", xout))]^2, na.rm=TRUE) / mean(mr_comp[, get(paste0("b_se_", yout))]^2, na.rm=TRUE)
# fit         <- mcr::mcreg(mr_comp[, get(paste0("b_", xout))], mr_comp[, get(paste0("b_", yout))], method.reg = "Deming", error.ratio = error_ratio, na.rm = TRUE)
# slope       <- coef(fit)[2]
# intercept   <- coef(fit)[1]
#
# ggplot(mr_comp, aes(x = get(paste0("b_", xout)), y = get(paste0("b_", yout)))) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_abline(intercept = intercept, slope = slope, color = "red") +
#   geom_text_repel(aes(label = exposure)) +
#   geom_point()
#


# p <- ggplot(hf, aes(x = b, y = exp_lab, color = outcome)) +
#   geom_point(position = position_dodge(width=0.5)) +
#   geom_errorbarh(aes(xmin = b-1.96*b_se, xmax=b+1.96*b_se), position = position_dodge(width=0.5)) +
#   scale_y_discrete(labels = function(x) sub("Pheno[0-9]_", "", x)) +
#   facet_wrap(~outcome, scales = "free_y", ncol=1)
# p

