#### DESCRIPTION ###############################################################
# observational HR vs genome-wide MR OR.

#### SET INPUT #################################################################
gw_mr_file    <- snakemake@input[["gw_mr"]]
name_map_file <- snakemake@input[["name_map"]]
liu_file      <- snakemake@input[["liu_fp"]]
julkunen_file <- snakemake@input[["julkunen_fp"]]
rep_file      <- snakemake@input[["rep"]]
out_plot      <- snakemake@output[["comparison_plot"]]
################################################################################

if (FALSE) {
  library(data.table); library(readxl)
  repo_dir      <- Sys.getenv("HF_METABOLITE_REPO2")
  gw_mr_file    <- file.path(repo_dir, "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz")
  name_map_file <- file.path(repo_dir, "scripts/gwas_metab_name_map.xlsx")
  liu_file      <- file.path(repo_dir, "scripts/circhf_circhf-2023-010896_supp2.xlsx")
  julkunen_file <- file.path(repo_dir, "scripts/julkunen_2023_summary_statistics.csv")
  rep_file      <- file.path(repo_dir, "output/tables/replication/metabolite_replication_bbs_direct_mr.tsv")
  out_plot      <- file.path(repo_dir, "output/figures/outcomes/observational_comparison.png")
}

message("[DEBUG] script start")
library(data.table)
library(ggplot2)
library(readxl)

dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

# setup
outcome_order  <- c("heart_failure", "hfref", "hfpef")
outcome_labels <- c(heart_failure = "All HF", hfref = "HFrEF", hfpef = "HFpEF")

# map
map_raw <- as.data.table(read_xlsx(name_map_file, sheet = 1))
liu_map <- map_raw[!is.na(liu_name)      & !is.na(gwas_id),
                   .(liu_name, gwas_id, label, pathway_group)]
jul_map <- map_raw[!is.na(julkunen_name) & !is.na(gwas_id),
                   .(julkunen_name, gwas_id, label, pathway_group)]

# bmi metabs
rep_dat  <- fread(rep_file)
bmi_labs <- rep_dat[replicates_trials_mr == "Trial & MR consistent", unique(label)]

# helper
parse_liu <- function(hr, ci_str) {
  hr_num <- suppressWarnings(as.numeric(hr))
  if (is.na(hr_num) || hr_num <= 0) return(list(b = NA_real_, b_se = NA_real_))
  parts <- suppressWarnings(
    as.numeric(trimws(strsplit(as.character(ci_str), ",")[[1]]))
  )
  if (length(parts) != 2 || any(is.na(parts)) || any(parts <= 0))
    return(list(b = log(hr_num), b_se = NA_real_))
  list(b    = log(hr_num),
       b_se = (log(parts[2]) - log(parts[1])) / (2 * 1.96))
}

# Liu et al.
hf_types <- c("all", "hfref", "hfpef")
col_vec  <- as.vector(unlist(lapply(hf_types, function(hf) {
  sapply(1:3, function(m) paste0(hf, "_model_", m, "_", c("hr", "ci_95", "p_fdr")))
})))
hf1 <- as.data.table(read_xlsx(liu_file, sheet = 1, skip = 4,
  col_names = c("liu_name", "pathway", "subpathway", "hmdb",
                col_vec, "new", "fig_num", "lasso")))
hf2 <- as.data.table(read_xlsx(liu_file, sheet = 2, skip = 4,
  col_names = c("liu_name", col_vec)))
liu <- rbind(hf1, hf2, fill = TRUE)[!is.na(liu_name)]

liu_long <- rbindlist(lapply(hf_types, function(hf_t) {
  sub <- liu[, .(
    liu_name,
    hr    = get(paste0(hf_t, "_model_1_hr")),
    ci_95 = get(paste0(hf_t, "_model_1_ci_95")),
    p_fdr = suppressWarnings(as.numeric(get(paste0(hf_t, "_model_1_p_fdr")))),
    outcome = fcase(hf_t == "all",   "heart_failure",
                    hf_t == "hfref", "hfref",
                    hf_t == "hfpef", "hfpef")
  )]
  parsed <- sub[, parse_liu(hr, ci_95), by = seq_len(nrow(sub))]
  sub[, .(liu_name, outcome, liu_b = parsed$b, liu_b_se = parsed$b_se, p_fdr)]
}))
liu_long <- liu_long[!is.na(liu_b) & !is.na(liu_b_se) & is.finite(liu_b) & is.finite(liu_b_se)]
liu_long <- merge(liu_long, liu_map, by = "liu_name", all.x = FALSE)
message("[DEBUG] liu_long rows: ", nrow(liu_long))

# Julkunen et al.
jul_raw  <- fread(julkunen_file)
jul_raw[, jul_b    := log(hazard_ratio)]
jul_raw[, jul_b_se := (log(ci_upper) - log(ci_lower)) / (2 * 1.96)]
jul_raw[, jul_p_fdr := p.adjust(pvalue, method = "fdr")]
jul_long <- jul_raw[is.finite(jul_b) & is.finite(jul_b_se),
                    .(julkunen_name = biomarker_name, jul_b, jul_b_se, jul_p_fdr)]
jul_long <- merge(jul_long, jul_map, by = "julkunen_name", all.x = FALSE)
message("[DEBUG] jul_long rows: ", nrow(jul_long))

# GW-MR
gw        <- fread(gw_mr_file)
gw        <- gw[method == "mr_ivw" & (is.na(error) | error == "") &
                outcome %in% outcome_order]
match_col <- if ("metabolite_id" %in% names(gw)) "metabolite_id" else "exposure"
gw_long   <- gw[, .(gwas_id = get(match_col), outcome, gw_b = b, gw_b_se = b_se, gw_p = p)]
message("[DEBUG] gw_long rows: ", nrow(gw_long))

# merge Liu
dat_liu <- merge(
  liu_long[label %in% bmi_labs],
  gw_long, by = c("gwas_id", "outcome"), all = FALSE
)
message("[DEBUG] dat_liu rows: ", nrow(dat_liu))

# merge Julkunen
dat_jul <- merge(
  jul_long[label %in% bmi_labs],
  gw_long[outcome == "heart_failure"], by = "gwas_id", all = FALSE
)
dat_jul[, outcome := "heart_failure"]
message("[DEBUG] dat_jul rows: ", nrow(dat_jul))

if (nrow(dat_liu) + nrow(dat_jul) == 0) {
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No overlapping data", hjust = 0.5) +
    theme_void()
  ggsave(out_plot, p_empty, width = 8, height = 4, dpi = 150, bg = "white")
  quit(save = "no", status = 0)
}

# OR cols
fmt_label <- function(x) {
  x <- trimws(gsub("_", " ", x))
  ifelse(nchar(x) <= 4, toupper(x), tools::toTitleCase(x))
}
add_cols <- function(d, b_col, b_se_col) {
  d[, obs_or   := exp(get(b_col))]
  d[, obs_lo95 := exp(get(b_col) - 1.96 * get(b_se_col))]
  d[, obs_hi95 := exp(get(b_col) + 1.96 * get(b_se_col))]
  d[, gw_or    := exp(gw_b)]
  d[, gw_lo95  := exp(gw_b - 1.96 * gw_b_se)]
  d[, gw_hi95  := exp(gw_b + 1.96 * gw_b_se)]
  d[, gw_sig   := !is.na(gw_p) & gw_p < 0.05]
  d[, metab_label := fmt_label(label)]
  d
}
dat_liu <- add_cols(dat_liu, "liu_b", "liu_b_se")
dat_jul <- add_cols(dat_jul, "jul_b", "jul_b_se")
dat_liu[, obs_sig := !is.na(p_fdr)     & p_fdr     < 0.05]
dat_jul[, obs_sig := !is.na(jul_p_fdr) & jul_p_fdr < 0.05]

# dont need X metabs
dat_liu <- dat_liu[!grepl("^X-\\d", label)]
dat_jul <- dat_jul[!grepl("^X-\\d", label)]
message("[DEBUG] after unknown filter — liu: ", nrow(dat_liu), "  jul: ", nrow(dat_jul))

# factors
all_pg <- sort(unique(c(dat_liu$pathway_group, dat_jul$pathway_group)))
dat_liu[, pathway_group := factor(pathway_group, levels = all_pg)]
dat_jul[, pathway_group := factor(pathway_group, levels = all_pg)]
message("[DEBUG] pathway groups: ", paste(all_pg, collapse = ", "))

# colours
n_pg    <- length(all_pg)
pg_cols <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(min(n_pg, 8), "Dark2"))(n_pg),
  all_pg
)

# combine
ci_clip_lo <- 0.15
ci_clip_hi <- 8.0
clip <- function(lo, hi) list(lo = pmax(lo, ci_clip_lo), hi = pmin(hi, ci_clip_hi))

gw_meta <- rbind(
  dat_liu[, .(gwas_id, outcome, metab_label, pathway_group, gw_sig,
              or   = gw_or,
              lo95 = clip(gw_lo95, gw_hi95)$lo,
              hi95 = clip(gw_lo95, gw_hi95)$hi)],
  dat_jul[!gwas_id %in% dat_liu[outcome == "heart_failure", gwas_id],
          .(gwas_id, outcome, metab_label, pathway_group, gw_sig,
            or   = gw_or,
            lo95 = clip(gw_lo95, gw_hi95)$lo,
            hi95 = clip(gw_lo95, gw_hi95)$hi)]
)
gw_meta <- unique(gw_meta, by = c("gwas_id", "outcome"))

long_dat <- rbind(
  dat_liu[, .(gwas_id, outcome, metab_label, pathway_group, gw_sig, obs_sig,
              or   = obs_or,
              lo95 = clip(obs_lo95, obs_hi95)$lo,
              hi95 = clip(obs_lo95, obs_hi95)$hi,
              source = "Liu et al. (obs)")],
  dat_jul[, .(gwas_id, outcome, metab_label, pathway_group, gw_sig, obs_sig,
              or   = obs_or,
              lo95 = clip(obs_lo95, obs_hi95)$lo,
              hi95 = clip(obs_lo95, obs_hi95)$hi,
              source = "Julkunen et al. (obs)")],
  gw_meta[, .(gwas_id, outcome, metab_label, pathway_group, gw_sig,
              obs_sig = FALSE,
              or, lo95, hi95,
              source = "GW-MR (IVW)")]
)
long_dat <- long_dat[is.finite(or)]
long_dat[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]

src_levels <- c("Liu et al. (obs)", "Julkunen et al. (obs)", "GW-MR (IVW)")
long_dat[, source := factor(source, levels = src_levels)]

# plotting types
pt_levels <- c("MS observational: FDR < 0.05", "MS observational: n.s.",
               "NMR observational: FDR < 0.05", "NMR observational: n.s.",
               "GW-MR: p < 0.05", "GW-MR: n.s.")
long_dat[, point_type := factor(fcase(
  source == "Liu et al. (obs)"      &  obs_sig, "MS observational: FDR < 0.05",
  source == "Liu et al. (obs)"      & !obs_sig, "MS observational: n.s.",
  source == "Julkunen et al. (obs)" &  obs_sig, "NMR observational: FDR < 0.05",
  source == "Julkunen et al. (obs)" & !obs_sig, "NMR observational: n.s.",
  source == "GW-MR (IVW)"          &  gw_sig,  "GW-MR: p < 0.05",
  source == "GW-MR (IVW)"          & !gw_sig,  "GW-MR: n.s."
), levels = pt_levels)]

# ordering
ord_hf     <- dat_liu[outcome == "heart_failure",
                      .(gwas_id, metab_label, pathway_group, sort_b = liu_b)]
ord_jul_hf <- dat_jul[
  !gwas_id %in% dat_liu[outcome == "heart_failure", gwas_id],
  .(gwas_id, metab_label, pathway_group, sort_b = jul_b)
]
ord_global <- unique(rbind(ord_hf, ord_jul_hf), by = "gwas_id")
metab_ord  <- ord_global[order(pathway_group, sort_b), unique(metab_label)]

long_dat[, metab_label := factor(metab_label, levels = metab_ord)]

# plot
p <- ggplot(long_dat, aes(x = or, y = metab_label, colour = pathway_group)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.4) +
  geom_line(aes(group = gwas_id), colour = "grey65", linewidth = 0.25, linetype = "dotted") +
  geom_errorbarh(
    aes(xmin = lo95, xmax = hi95),
    height = 0
  ) +
  geom_point(aes(shape = point_type, size = point_type,
                 fill = pathway_group), stroke = 0.4, alpha = 0.9) +
  scale_shape_manual(name = "Estimate",
    values = c("MS observational: FDR < 0.05"  = 25,      # filled down triangle
               "MS observational: n.s."         = 6,   # open down triangle
               "NMR observational: FDR < 0.05" = 24,  # filled up triangle
               "NMR observational: n.s."        = 2,   # open up triangle
               "GW-MR: p < 0.05"               = 16,  # filled circle
               "GW-MR: n.s."                   = 1))  + # open circle
  scale_size_manual(name = "Estimate",
    values = c("MS observational: FDR < 0.05"  = 1.5,
               "MS observational: n.s."         = 1.5,
               "NMR observational: FDR < 0.05" = 1.5,
               "NMR observational: n.s."        = 1.5,
               "GW-MR: p < 0.05"               = 2.0,
               "GW-MR: n.s."                   = 2.0)) +
  scale_colour_manual(name = "Pathway", values = pg_cols) +
  scale_fill_manual(values = pg_cols, guide = "none") +
  guides(shape = guide_legend(
    title.position = "top",
    override.aes   = list(
      fill   = c("grey30", NA, "grey30", NA, "grey30", NA),
      colour = rep("grey30", 6),
      size   = 3
    )
  )) +
  scale_x_log10(
    breaks = function(x) {
      cands <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
      cands[cands >= x[1] & cands <= x[2]]
    },
    labels = function(x) ifelse(x == round(x), as.integer(x), x)
  ) +
  facet_wrap(~ outcome_label, nrow = 1, scales = "free_x") +
  labs(
    x = "HR / OR (95% CI, log scale)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey70"),
    legend.position  = "right",
    legend.key.size  = unit(0.9, "lines"),
    legend.text      = element_text(size = 8.5),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    plot.caption     = element_text(size = 7.5, hjust = 0, colour = "grey40"),
    plot.margin      = margin(8, 8, 8, 8),
    axis.text.y      = element_text(size = 8)
  )


ggsave(out_plot, p,
       width  = 9.5,
       height = 10,
       dpi    = 300, bg = "white")
message("[DEBUG] saved: ", out_plot)
