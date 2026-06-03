#### DESCRIPTION #######################################
# Heatmap of metabolite → outcome MR results
########################################################

#### SET INPUT #########################################
results_file     <- snakemake@input[["results_file"]]
rev_mr_file      <- snakemake@input[["rev_mr_file"]]
name_map_file    <- snakemake@input[["name_map_file"]]
replicating_file <- snakemake@input[["replicating_file"]]
cluster_file     <- snakemake@input[["cluster_file"]]
output_file      <- snakemake@output[["heatmap"]]
log_file         <- snakemake@log[["log"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
# fl <- list.files(file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tmp_objects/mr_results"), pattern = "_on_(hip_knee_osteoarthritis|heart_failure|endometrial_cancer)_mr\\.tsv", full.names = T)
# d  <- rbindlist(lapply(fl, fread), fill = T)
results_file <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/metabolites_on_outcomes.tsv")
# fwrite(d, results_file, sep="\t")
# fl_rev <- list.files(file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tmp_objects/mr_results"), pattern = "hip_knee_osteoarthritis|heart_failure|endometrial_cancer)_on_GCST[0-9]+_mr\\.tsv", full.names = T)
# d_rev  <- rbindlist(lapply(fl_rev, fread), fill = T)
rev_mr_file <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/outcomes_on_metabolites.tsv")
# fwrite(d_rev, rev_mr_file, sep="\t")
name_map_file    <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "scripts/gwas_metab_name_map.xlsx")
replicating_file <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/replication/study_replication.tsv")
cluster_file     <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/instruments/cluster_membership.tsv")
output_file      <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/figures/outcomes/metab_outcome_mr.png")
log_file         <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/logs/metab_outcome_mr_instruments.log")
formatted_mr_results <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/formatted_metabolites_on_outcomes.tsv")
# ─────────────────────────────────────────────────────────────────────────────
}

library(data.table)
library(ggplot2)
library(readxl)

dt     <- fread(results_file, sep = "\t")
rev_mr <- fread(rev_mr_file)
clust  <- fread(cluster_file, sep = "\t")
rep    <- fread(replicating_file)
map    <- read_xlsx(name_map_file, sheet=1) |> as.data.table()

# disambiguate labels that appear on more than one platform
platform_suffix <- c("nightingale" = "NMR", "metabolon" = "MS")
dup_labels      <- map[, .(n = uniqueN(gwas_id)), by = label][n > 1, label]
map[label %in% dup_labels,
    label := paste0(label, " (", platform_suffix[platform], ")")]
rep[label %in% dup_labels,
    label := paste0(label, " (", platform_suffix[platform], ")")]

# annotate
dt[map, c("label","pathway_group","platform") := .(i.label, i.pathway_group, i.platform), on = c("exp_id"="gwas_id")]
dt[, label := fcoalesce(label, exp_id)]
dt[rep, replicates := i.replicates, on = "label"]
dt <- dt[replicates == "Trial consistent (interaction)"]
dt[clust, c("cluster_id","cluster_singleton") := .(i.cluster_id, i.cluster_singleton), on = c("exp_id"="gwas_id")]
dt <- dt[!is.na(cluster_id),]
clust[unique(dt[, .(gwas_id = exp_id, label)]), label := i.label, on = "gwas_id"]

# adjust pvalue
dt[, p_fdr := p.adjust(p, method = "fdr", n = .N), by = .(out_id, method, steiger_filtering)]
setcolorder(dt, c("p_fdr"), after = "p")

# write out forward MR
exp_meta <- dt[steiger_filtering == FALSE, .(
  label      = label[1],
  pathway_group = pathway_group[1],
  platform   = platform[1],
  cluster_id = cluster_id[1],
  cluster_singleton = cluster_singleton[1],
  fstat      = fstat[1],
  exp_snp_r2 = exp_snp_r2[1],
  out_snp_r2 = out_snp_r2[1]
), by = exp_id]

dt_formatted <- dt[steiger_filtering == FALSE, .(
  exp_id, out_id, method,
  nsnp,
  b, b_se,
  p, p_fdr, Q, Q_pval,
  egger_intercept, egger_intercept_se, egger_intercept_p,
  presso_error, presso_global_p, n_outliers, prop_outliers, distortion_p,
  steiger_correct_dir, steiger_p
)]

complete_grid <- CJ(
  exp_id = dt_formatted[, unique(exp_id)],
  out_id = dt_formatted[, unique(out_id)],
  method = dt_formatted[, unique(method)]
)

dt_formatted <- dt_formatted[complete_grid, on = .(exp_id, out_id, method)]
dt_formatted <- exp_meta[dt_formatted, on = "exp_id"]
dt_formatted[is.na(b) | !grepl("Inverse variance weighted|Wald ratio", method),
             c("steiger_correct_dir","steiger_p","exp_snp_r2","out_snp_r2") := NA]
dt_formatted[, or := ifelse(is.na(b),
                            NA_character_,
                            sprintf("%.3f (%.3f-%.3f)", exp(b), exp(b - 1.96*b_se), exp(b + 1.96*b_se)))]
setcolorder(dt_formatted, "or", before = "b")
dt_formatted[, c("b","b_se") := NULL]
setcolorder(dt_formatted, c("exp_snp_r2","out_snp_r2"), after = "steiger_p")
dt_formatted[, method := factor(method,
  levels = c("Wald ratio", "Inverse variance weighted", "MR Egger",
             "Weighted median", "Weighted mode", "MR-PRESSO-raw", "MR-PRESSO-corrected")
)]
dt_formatted <- dt_formatted[order(out_id, exp_id, method)]
#dt_formatted <- dt_formatted[, lapply(.SD, function(x) ifelse(is.na(x) | x=="", "-", as.character(x)))]

fwrite(dt_formatted, formatted_mr_results, sep="\t")



# Rename Wald ratio → IVW and append Steiger-filtered IVW as a separate method level
merge_steiger <- function(d) {
  d[method == "Wald ratio", method := "Inverse variance weighted"]
  st <- d[steiger_filtering == TRUE & method == "Inverse variance weighted"]
  st[, method := "Inverse variance weighted Steiger filtered"]
  rbind(d[steiger_filtering == FALSE], st)[order(exp_id, out_id, method)]
}

# 6-method per-outcome stream — built before dt is reduced to 5 methods below
dt_per_out <- merge_steiger(copy(dt))
dt_per_out[, out_label := tools::toTitleCase(gsub("_", " ", out_id))]

dt     <- merge_steiger(dt)
rev_mr <- merge_steiger(rev_mr)

method_levels <- c(
  "Inverse variance weighted",
  "MR-PRESSO-corrected",
  "Weighted median",
  "Weighted mode",
  "MR Egger"
)
method_labels <- c(
  "Inverse variance weighted" = "IVW",
  "MR-PRESSO-corrected"       = "MR-PRESSO",
  "Weighted median"           = "Weighted median",
  "Weighted mode"             = "Weighted mode",
  "MR Egger"                  = "MR-Egger"
)

dt <- dt[method %in% method_levels]
dt[, method := factor(method, levels = method_levels)]

# signed -log10(p): positive = increases risk, negative = reduces risk
dt[, signed_log10p := sign(b) * (-log10(p))]

# outcome labels
dt[, out_label := tools::toTitleCase(gsub("_", " ", out_id))]

# per-label mean IVW -log10(p) for ordering within each cluster
ivw_ord <- dt[method == "Inverse variance weighted",
              .(mean_log10p = mean(-log10(p), na.rm = TRUE)),
              by = label]

ivw_ord[clust, cluster := i.cluster, on = "label"]
ivw_ord[is.na(cluster), cluster := 0L]

non_sing_clusters <- sort(unique(ivw_ord[cluster > 0L, cluster]))
label_levels <- c(
  unlist(lapply(non_sing_clusters, function(cl)
    ivw_ord[cluster == cl][order(mean_log10p), as.character(label)]
  )),
  ivw_ord[cluster == 0L][order(mean_log10p), as.character(label)]
)
dt[, label := factor(label, levels = label_levels)]

# cluster_group column for faceting
cluster_group_levels <- c(paste0("Cluster ", non_sing_clusters), "Singletons")
dt[ivw_ord, cluster := i.cluster, on = "label"]
dt[, cluster_group := factor(
  fifelse(cluster == 0L, "Singletons", paste0("Cluster ", cluster)),
  levels = cluster_group_levels
)]

# sector colours for inner circos track — derived from the cluster colour table
sector_col_dt <- clust[, .(
  cluster_group = fifelse(cluster == 0L, "Singletons", paste0("Cluster ", cluster)),
  color         = fifelse(cluster == 0L, "grey80", color)
)][, .(color = color[1L]), by = cluster_group]
sector_cols <- setNames(sector_col_dt$color, sector_col_dt$cluster_group)[cluster_group_levels]

# instrument info panel: nsnp + fstat per label × out_label, from IVW row
dt_info <- dt[method == "Inverse variance weighted",
              .(nsnp  = nsnp[1L],
                fstat = if ("fstat" %in% names(dt)) fstat[1L] else NA_real_),
              by = .(label, out_label, cluster_group)]
dt_info[, method    := "Instruments"]
dt_info[, info_label := paste0(nsnp, "/", ifelse(is.na(fstat), "?", round(fstat, 0)))]

# extend factor levels and combine — Instruments first
all_levels <- c("Instruments", levels(dt$method))
dt[,      method := factor(as.character(method), levels = all_levels)]
dt_info[, method := factor(method,               levels = all_levels)]
dt_info[, label  := factor(label,                levels = levels(dt$label))]
dt_plot <- rbind(dt, dt_info, fill = TRUE)
dt_plot[, method := factor(as.character(method), levels = all_levels)]

method_labels["Instruments"] <- "N SNPs / F"

# ── Instrument summary log ────────────────────────────────────────────────────
{
  rep[map, gwas_id := i.gwas_id, on = "label"]
  n_rep_total        <- uniqueN(rep[replicates == "Trial consistent (interaction)", .(label, platform)])
  n_with_gwas_ids    <- rep[replicates == "Trial consistent (interaction)" & !is.na(gwas_id), .N]
  n_with_instruments <- dt[method == "Inverse variance weighted", uniqueN(label)]
  pct_instruments    <- round(100 * n_with_instruments / n_rep_total, 1)
  n_with_instruments_post_harm <- dt[method == "Inverse variance weighted" & nsnp > 0, uniqueN(label)]
  pct_instruments_post_harm    <- round(100 * n_with_instruments_post_harm / n_rep_total, 1)


  nsnp_vals  <- dt_info[!is.na(nsnp),  nsnp]
  fstat_vals <- dt_info[!is.na(fstat), fstat]

  fmt_stats <- function(vals) {
    sprintf(
      "range %d-%d, median %d, IQR %d-%d",
      as.integer(min(vals)), as.integer(max(vals)),
      as.integer(median(vals)),
      as.integer(quantile(vals, 0.25)), as.integer(quantile(vals, 0.75))
    )
  }

  ivw_dt <- dt[as.character(method) == "Inverse variance weighted"]
  out_sig <- ivw_dt[, .(
    n_total  = uniqueN(label),
    n_nom    = sum(p     < 0.05, na.rm = TRUE),
    n_fdr    = sum(p_fdr < 0.05, na.rm = TRUE)
  ), by = out_label][order(out_label)]

  out_sig_lines <- apply(out_sig, 1, function(r) {
    sprintf("  %-35s total=%s  p<0.05: %s (%.1f%%)  FDR p<0.05: %s (%.1f%%)",
            r["out_label"],
            r["n_total"],
            r["n_nom"], 100 * as.integer(r["n_nom"]) / as.integer(r["n_total"]),
            r["n_fdr"], 100 * as.integer(r["n_fdr"]) / as.integer(r["n_total"]))
  })

  log_lines <- c(
    "=== Metabolite -> Outcome MR: Instrument Summary ===",
    "",
    sprintf("a) Replicating metabolites (Trial consistent) in replication file:  %d", n_rep_total),
    sprintf("b) Replicating metabolites (Trial consistent) with a GWAS ID:       %d", n_with_gwas_ids),
    sprintf("c) Of these with genetic instruments:                               %d (%.1f%%)", n_with_instruments, pct_instruments),
    sprintf("d) Of these with valid genetic instruments (in MR analysis):        %d (%.1f%%)", n_with_instruments_post_harm, pct_instruments_post_harm),
    "",
    "--- IVW nSNP distribution (across all metabolite x outcome analyses) ---",
    sprintf("    %s", fmt_stats(nsnp_vals)),
    "",
    "--- IVW F-statistic distribution (across all metabolite x outcome analyses) ---",
    sprintf("    %s", fmt_stats(fstat_vals)),
    "",
    "--- IVW significant associations per outcome ---",
    out_sig_lines
  )

  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(log_lines, log_file)
}

# cap colour scale at 99th percentile to prevent outliers dominating
abs_max <- quantile(abs(dt$signed_log10p), 0.99, na.rm = TRUE)

p <- ggplot(dt_plot, aes(x = method, y = label)) +
  geom_tile(
    data      = dt_plot[as.character(method) != "Instruments"],
    aes(fill  = signed_log10p),
    colour    = "white", linewidth = 0.3, na.rm = TRUE
  ) +
  geom_tile(
    data   = dt_info,
    fill   = "grey95", colour = "white", linewidth = 0.3
  ) +
  geom_text(
    data   = dt_plot[!is.na(p_fdr) & p_fdr < 0.05],
    aes(label = "x"),
    colour = "grey10", size = 2.6, vjust = 0.5, na.rm = TRUE
  ) +
  geom_text(
    data  = dt_info,
    aes(label = info_label),
    size  = 1.75, colour = "grey30"
  ) +
  scale_fill_gradient2(
    name     = expression(paste("Signed ", -log[10](italic(p)))),
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-abs_max, abs_max),
    oob      = scales::squish,
    na.value = "grey88"
  ) +
  scale_x_discrete(limits = all_levels, labels = method_labels) +
  facet_grid(cluster_group ~ out_label, scales = "free", space = "free_y") +
  labs(x = NULL, y = "Metabolite (exposure)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x       = element_text(angle = 40, hjust = 1, size = 9),
    axis.text.y       = element_text(size = 8),
    axis.ticks        = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    panel.spacing.x   = unit(0.5, "lines"),
    strip.background  = element_rect(fill = "grey20", colour = NA),
    strip.text        = element_text(colour = "white", face = "bold", size = 9),
    legend.position   = "right",
    legend.key.height = unit(1.8, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.title      = element_text(size = 9),
    plot.margin       = margin(8, 8, 8, 8)
  )
p

n_metabs   <- uniqueN(dt$label)
n_methods  <- uniqueN(dt$method)   # real MR methods only
n_outcomes <- uniqueN(dt$out_label)

fig_w <- max(3, 1.0 + ((n_methods + 1) * 0.45 + 0.3) * n_outcomes)
fig_h <- max(2, 0.6 + n_metabs * 0.15)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, p, width = fig_w, height = fig_h, dpi = 300, bg = "white", units = "in")
cat("Saved:", output_file, "\n")






# ══ CIRCOS PLOT ══════════════════════════════════════════════════════════════
library(circlize)

metab_order <- levels(dt$label)
lookup_dt   <- unique(dt[, .(label, cluster_group)])
pathway_split <- factor(
  as.character(lookup_dt[match(metab_order, as.character(lookup_dt$label)), cluster_group]),
  levels = cluster_group_levels
)
length(pathway_split)
length(metab_order)
table(pathway_split, useNA = "ifany")

gap <- 18
gap_after_vals <- setNames(rep(3, nlevels(pathway_split)), levels(pathway_split))
gap_after_vals[levels(pathway_split)[nlevels(pathway_split)]] <- gap

# explicit outer→inner order
outcomes_vec  <- intersect(
  c("Heart Failure", "Hip Knee Osteoarthritis", "Kidney Cancer" ),
  unique(dt$out_label)
)
n_out         <- length(outcomes_vec)
n_meth        <- length(method_levels)
last_sector   <- levels(pathway_split)[length(levels(pathway_split))]

# Discrete colour scheme: direction × significance tier (IVW FDR as primary test)
col_blue <- "#2166AC"
col_red  <- "#B2182B"
c_blue_fdr <- adjustcolor(col_blue, alpha.f = 1.00)
c_blue_nom <- adjustcolor(col_blue, alpha.f = 0.50)
c_blue_ns  <- adjustcolor(col_blue, alpha.f = 0.15)
c_red_fdr  <- adjustcolor(col_red,  alpha.f = 1.00)
c_red_nom  <- adjustcolor(col_red,  alpha.f = 0.50)
c_red_ns   <- adjustcolor(col_red,  alpha.f = 0.15)

sig_colour <- function(b, p, p_fdr) {
  fcase(
    is.na(p) | is.na(b),   "white",
    b < 0 & p_fdr < 0.05,  c_blue_fdr,
    b < 0 & p < 0.05,       c_blue_nom,
    b < 0,                   c_blue_ns,
    b > 0 & p_fdr < 0.05,  c_red_fdr,
    b > 0 & p < 0.05,       c_red_nom,
    b > 0,                   c_red_ns,
    default                 = "white"
  )
}

# One character colour matrix per outcome: rows = metabolites, cols = MR methods.
# sig_colour is computed in long format so b, p, p_fdr_ivw always come from the
# same row — avoids any column-ordering mismatch from separate dcasts.
make_circ_mat <- function(out) {
  sub_dt <- dt[out_label == out & as.character(method) %in% method_levels,
               .(label, method, b, p, p_fdr)]
  sub_dt[, col := sig_colour(b, p, p_fdr)]
  wide   <- dcast(sub_dt, label ~ method, value.var = "col", fun.aggregate = function(x) x[1L])
  ord    <- match(metab_order, as.character(wide$label))
  m_col  <- as.matrix(wide[ord, method_levels, with = FALSE])
  m_col[is.na(m_col)] <- "white"
  rownames(m_col) <- metab_order
  colnames(m_col) <- unname(method_labels[method_levels])
  m_col
}
circ_mats <- setNames(lapply(outcomes_vec, make_circ_mat), outcomes_vec)

all_col_vals <- unique(unlist(lapply(circ_mats, as.vector)))
col_disc     <- setNames(all_col_vals, all_col_vals)

# Per-metabolite average nSNP (across outcomes that actually had an instrument) and Fstat across outcomes (for instrument track)
metab_nsnp  <- dt_info[, .(avg_nsnp  = mean(nsnp[!is.na(fstat)],  na.rm = TRUE)), by = label]
metab_fstat <- dt_info[, .(avg_fstat = mean(fstat,  na.rm = TRUE)), by = label]

nsnp_vec  <- metab_nsnp[ match(metab_order, as.character(label)), avg_nsnp]
fstat_vec <- metab_fstat[match(metab_order, as.character(label)), avg_fstat]

sector_mets <- split(metab_order, pathway_split)

# Egger intercept p < 0.05 indicator: named list per outcome, TRUE/FALSE per metab
egger_intercept_sig <- lapply(outcomes_vec, function(out) {
  sub <- dt[out_label == out & method == "MR Egger",
            .(label = as.character(label), egger_intercept_p)]
  setNames(
    !is.na(sub$egger_intercept_p) & sub$egger_intercept_p < 0.05,
    sub$label
  )
})
names(egger_intercept_sig) <- outcomes_vec

# y-centre of the MR Egger column within each heatmap track
# column j displayed at y = n_meth - j + 0.5; Egger is column 4
egger_y <- n_meth - match("MR Egger", method_levels) + 0.5

circos_file <- sub("\\.png$", "_circos.png", output_file)
meth_abbr   <- rev(c("IVW", "PRESSO", "W.Med", "W.Mode", "Egger"))

png(circos_file, width = 10, height = 10, units = "in", res = 600, bg = "white")
circos.clear()
circos.par(
  gap.after               = gap_after_vals,
  track.margin            = c(0.005, 0.018),  # 0.018 outer margin = space for label band
  cell.padding            = c(0, 0, 0, 0),
  start.degree            = 270 - gap,
  points.overflow.warning = FALSE,
  canvas.xlim             = c(-1.05, 1.05),
  canvas.ylim             = c(-1.05, 1.05)
)

# Cluster membership track — outermost ring, initialises the circle.
# One thin coloured band per metabolite row using the same palette as the
# instrument-overlap graph.
clust_col_vec           <- setNames(clust$color, clust$label)
clust_char_vec          <- clust_col_vec[metab_order]
clust_char_vec[is.na(clust_char_vec)] <- "grey80"
clust_mat               <- matrix(clust_char_vec, ncol = 1,
                                  dimnames = list(metab_order, ""))
uniq_cols               <- unique(clust_char_vec)
clust_col_map           <- setNames(uniq_cols, uniq_cols)

circos.heatmap(
  clust_mat,
  split              = pathway_split,
  col                = clust_col_map,
  cluster            = FALSE,
  track.height       = 0.015,
  track.margin       = c(0.005, 0.005),
  show.sector.labels = FALSE,
  rownames.side      = "outside",
  rownames.cex       = 0.3,
  rownames.font      = 1,
  cell.border        = NA
)
clust_tidx <- get.current.track.index()

for (sec in names(sector_mets)) {
  n_in_sec <- length(sector_mets[[sec]])
  if (n_in_sec == 0) next
  circos.text(
    x = n_in_sec / 2, y = 0.5,
    labels = if (sec == "Singletons") "Singletons" else paste0("C", sub("^Cluster ", "", sec)),
    sector.index = sec, track.index = clust_tidx,
    facing = "bending.inside", niceFacing = TRUE,
    cex = 0.25, col = "white", adj = c(0.5, 0.5)
  )
}

x_right_clust <- get.cell.meta.data("xlim", sector.index = last_sector,
                                     track.index = clust_tidx)[2]
circos.text(
  x            = x_right_clust,
  y            = 0.5,
  labels       = "Instrument Cluster ",
  sector.index = last_sector,
  track.index  = clust_tidx,
  facing       = "inside",
  niceFacing   = TRUE,
  cex          = 0.3,
  adj          = c(0, 0.5)
)

# Average nSNP track — text values
circos.track(
  ylim         = c(0, 1),
  track.height = 0.012,
  track.margin = c(0.002, 0.002),
  bg.border    = NA,
  bg.col       = NA,
  panel.fun    = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    mets <- sector_mets[[sec]]
    vals <- nsnp_vec[match(mets, metab_order)]
    for (k in seq_along(mets)) {
      if (!is.na(vals[k]))
        circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
    }
  }
)
nsnp_tidx    <- get.current.track.index()
x_right_nsnp <- get.cell.meta.data("xlim", sector.index = last_sector,
                                    track.index = nsnp_tidx)[2]
circos.text(
  x = x_right_nsnp, y = 0.5, labels = "Avg nSNP ",
  sector.index = last_sector, track.index = nsnp_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)

# Average Fstat track — text values
circos.track(
  ylim         = c(0, 1),
  track.height = 0.012,
  track.margin = c(0.002, 0.002),
  bg.border    = NA,
  bg.col       = NA,
  panel.fun    = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    mets <- sector_mets[[sec]]
    vals <- fstat_vec[match(mets, metab_order)]
    for (k in seq_along(mets)) {
      if (!is.na(vals[k]))
        circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
    }
  }
)
fstat_tidx    <- get.current.track.index()
x_right_fstat <- get.cell.meta.data("xlim", sector.index = last_sector,
                                     track.index = fstat_tidx)[2]
circos.text(
  x = x_right_fstat, y = 0.5, labels = "Avg Fstat ",
  sector.index = last_sector, track.index = fstat_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)

# Heatmap rings — circle already initialised above by the cluster track.
# Outcome labels drawn afterwards via draw.sector() in the outer track margin.
heatmap_tidx <- integer(length(outcomes_vec))
first_sec    <- levels(pathway_split)[1]
for (i in seq_along(outcomes_vec)) {
  circos.heatmap(
    circ_mats[[outcomes_vec[i]]],
    split              = pathway_split,
    col                = col_disc,
    na.col             = "white",
    cluster            = FALSE,
    cell.border        = "gray95",
    cell.lwd           = 0.5,
    track.height       = 0.11,
    show.sector.labels = FALSE,
    rownames.side      = "none",
    rownames.cex       = 0.3,
    rownames.font      = 1
  )
  heatmap_tidx[i] <- get.current.track.index()
}

# Overlay "*" on Egger cells where intercept p < 0.05 (directional pleiotropy signal)
for (i in seq_along(outcomes_vec)) {
  sig_map <- egger_intercept_sig[[outcomes_vec[i]]]
  for (sec in names(sector_mets)) {
    mets <- sector_mets[[sec]]
    for (k in seq_along(mets)) {
      if (isTRUE(sig_map[mets[k]])) {
        circos.text(
          x            = k - 0.5,
          y            = egger_y,
          labels       = "*",
          sector.index = sec,
          track.index  = heatmap_tidx[i],
          facing       = "inside",
          niceFacing   = TRUE,
          cex          = 0.4,
          col          = "grey40",
          adj          = c(0.5, 0.5)
        )
      }
    }
  }
}


# Method labels in the 15-degree gap after the last sector.
# facing="clockwise" keeps text within the radial bounds of each track
# (no bleed into neighbouring rings), extending into the angular gap.
x_right <- get.cell.meta.data("xlim", sector.index = last_sector,
                               track.index = heatmap_tidx[1])[2]
for (i in seq_along(outcomes_vec)) {
  for (j in seq_len(n_meth)) {
    circos.text(
      x            = x_right,
      y            = j - 0.5,
      labels       = paste0(meth_abbr[j], " "),
      sector.index = last_sector,
      track.index  = heatmap_tidx[i],
      facing       = "inside",
      niceFacing   = TRUE,
      cex          = 0.3,
      adj          = c(0, 0.5)
    )
  }
}

# Outcome labels — white arc in the outer track-margin of each heatmap ring,
# with a single text label at 12 o'clock (x=0, y=r in circlize canvas coords).
for (i in seq_along(outcomes_vec)) {
  top_r <- get.cell.meta.data("cell.top.radius",
                               sector.index = first_sec,
                               track.index  = heatmap_tidx[i])
  draw.sector(
    start.degree = 0, end.degree = 360,
    rou1 = top_r + 0.017, rou2 = top_r + 0.001,
    col = "white", border = NA
  )
  xlim1 <- get.cell.meta.data("xlim", sector.index = first_sec, track.index = heatmap_tidx[i])
  # y in data coords for the midpoint of the outer margin:
  # outer_margin_radial / track_height * n_meth = 0.009/0.13 * n_meth
  y_label <- n_meth + 0.013 / 0.13 * n_meth
  op <- par(xpd = NA)   # disable clipping so text outside ylim is still drawn
  circos.text(
    x            = 90 + 28,
    y            = y_label,
    labels       = outcomes_vec[i],
    sector.index = first_sec,
    track.index  = heatmap_tidx[i],
    facing       = "bending.inside",
    niceFacing   = TRUE,
    cex          = 0.4,
    font         = 2,
    col          = "grey20",
    adj          = c(0.5, 0.5)
  )
  par(op)
}

legend("topright", inset = c(0.1, 0.15),
       legend = c("FDR p<0.05 (risk)", "p<0.05 (risk)", "p≥0.05 (risk)",
                  "p≥0.05 (protective)", "p<0.05 (protective)", "FDR p<0.05 (protective)"),
       fill   = c(c_red_fdr, c_red_nom, c_red_ns, c_blue_ns, c_blue_nom, c_blue_fdr),
       title  = "MR result (FDR within outcome × method)", bty = "n", cex = 0.55)

circos.clear()
dev.off()
cat("Saved:", circos_file, "\n")






# ══ PER-OUTCOME CIRCOS PLOTS ══════════════════════════════════════════════════
# One plot per outcome; same tracks as the combined plot above, plus a reverse
# MR (outcome → metabolite) heatmap ring inside the forward MR ring.

# Label lookup: reuse disambiguated labels already applied to dt
label_lookup  <- unique(dt[, .(gwas_id = exp_id, label)])
out_id_lookup <- unique(dt[, .(out_label, out_id)])

# 6-method setup for per-outcome plots (IVW Steiger in both forward and reverse)
method_levels_per_out <- c(
  "Inverse variance weighted",
  "Inverse variance weighted Steiger filtered",
  "MR-PRESSO-corrected",
  "Weighted median",
  "Weighted mode",
  "MR Egger"
)
method_labels_per_out <- c(
  "Inverse variance weighted"                  = "IVW",
  "Inverse variance weighted Steiger filtered" = "IVW Steiger",
  "MR-PRESSO-corrected"                        = "MR-PRESSO",
  "Weighted median"                            = "Weighted median",
  "Weighted mode"                              = "Weighted mode",
  "MR Egger"                                   = "MR-Egger"
)
n_meth_per_out    <- length(method_levels_per_out)
meth_abbr_per_out <- rev(c("IVW", "IVW Steiger", "PRESSO", "W.Med", "W.Mode", "Egger"))
egger_y_per_out   <- n_meth_per_out - match("MR Egger", method_levels_per_out) + 0.5

rev_mr <- rev_mr[method %in% method_levels_per_out]
rev_mr[, method := factor(method, levels = method_levels_per_out)]
rev_mr[label_lookup, label := i.label, on = c("out_id" = "gwas_id")]
rev_mr[, p_fdr_rev := p.adjust(p, method = "fdr"), by = .(exp_id, method)]

dt_per_out <- dt_per_out[method %in% method_levels_per_out]
dt_per_out[, method := factor(method, levels = method_levels_per_out)]


for (this_outcome in outcomes_vec) {

  this_out_id     <- out_id_lookup[out_label == this_outcome, out_id][1]
  out_circos_file <- sub("\\.png$", paste0("_circos_", gsub("[^a-z0-9]", "_", tolower(this_outcome)), ".png"), output_file)

  # forward MR matrix — rebuilt from dt_per_out with 6 methods
  sub_fwd <- dt_per_out[out_label == this_outcome]
  sub_fwd[, col := sig_colour(b, p, p_fdr)]
  wide_fwd     <- dcast(sub_fwd, label ~ method, value.var = "col", fun.aggregate = function(x) x[1L])
  ord_fwd      <- match(metab_order, as.character(wide_fwd$label))
  circ_mat_fwd <- as.matrix(wide_fwd[ord_fwd, method_levels_per_out, with = FALSE])
  circ_mat_fwd[is.na(circ_mat_fwd)] <- "white"
  rownames(circ_mat_fwd) <- metab_order
  colnames(circ_mat_fwd) <- unname(method_labels_per_out[method_levels_per_out])

  # reverse MR matrix
  rev_mat  <- NULL
  rev_this <- rev_mr[exp_id == this_out_id & !is.na(label)]
  if (nrow(rev_this) > 0) {
    rev_this[, col := sig_colour(b, p, p_fdr_rev)]
    wide_rev <- dcast(rev_this, label ~ method, value.var = "col", fun.aggregate = function(x) x[1L])
    ord_rev  <- match(metab_order, as.character(wide_rev$label))
    rev_mat  <- as.matrix(wide_rev[ord_rev, method_levels_per_out, with = FALSE])
    rev_mat[is.na(rev_mat)] <- "white"
    rownames(rev_mat) <- metab_order
    colnames(rev_mat) <- unname(method_labels_per_out[method_levels_per_out])
  }

  # forward nSNP / Fstat / Steiger-excluded nSNP
  dt_info_fwd_i <- dt_per_out[method == "Inverse variance weighted" & out_label == this_outcome,
                               .(nsnp  = nsnp[1L],
                                 fstat = if ("fstat" %in% names(dt_per_out)) fstat[1L] else NA_real_),
                               by = label]
  nsnp_fwd_i  <- dt_info_fwd_i[match(metab_order, as.character(label)), nsnp]
  fstat_fwd_i <- dt_info_fwd_i[match(metab_order, as.character(label)), fstat]
  steiger_nsnp_fwd_dt <- dt_per_out[
    method == "Inverse variance weighted Steiger filtered" & out_label == this_outcome,
    .(steiger_nsnp = nsnp[1L]), by = label]
  steiger_nsnp_fwd_i <- steiger_nsnp_fwd_dt[match(metab_order, as.character(label)), steiger_nsnp]

  # reverse nSNP / Fstat / Steiger-excluded nSNP
  nsnp_rev_i         <- NULL
  fstat_rev_i        <- NULL
  steiger_nsnp_rev_i <- NULL
  if (!is.null(rev_mat) && nrow(rev_this) > 0) {
    rev_info_i <- rev_this[method == "Inverse variance weighted",
                           .(nsnp  = nsnp[1L],
                             fstat = if ("fstat" %in% names(rev_this)) fstat[1L] else NA_real_),
                           by = label]
    nsnp_rev_i  <- rev_info_i[match(metab_order, as.character(label)), nsnp]
    fstat_rev_i <- rev_info_i[match(metab_order, as.character(label)), fstat]
    rev_steiger_dt <- rev_this[method == "Inverse variance weighted Steiger filtered",
                               .(steiger_nsnp = nsnp[1L]), by = label]
    steiger_nsnp_rev_i <- rev_steiger_dt[match(metab_order, as.character(label)), steiger_nsnp]
  }

  # unified colour map for this plot
  all_cols_i <- unique(c(as.vector(circ_mat_fwd), if (!is.null(rev_mat)) as.vector(rev_mat), "white"))
  col_disc_i <- setNames(all_cols_i, all_cols_i)

  png(out_circos_file, width = 10, height = 10, units = "in", res = 600, bg = "white")
  circos.clear()
  circos.par(
    gap.after               = gap_after_vals,
    track.margin            = c(0.005, 0.018),
    cell.padding            = c(0, 0, 0, 0),
    start.degree            = 270 - gap,
    points.overflow.warning = FALSE,
    canvas.xlim             = c(-1.05, 1.05),
    canvas.ylim             = c(-1.05, 1.05)
  )

  # ── Cluster track ─────────────────────────────────────────────────────────
  circos.heatmap(
    clust_mat, split = pathway_split, col = clust_col_map,
    cluster = FALSE, track.height = 0.015, track.margin = c(0.005, 0.005),
    show.sector.labels = FALSE, rownames.side = "outside",
    rownames.cex = 0.3, rownames.font = 1, cell.border = NA
  )
  clust_ti <- get.current.track.index()

  for (sec in names(sector_mets)) {
    n_in_sec <- length(sector_mets[[sec]])
    if (n_in_sec == 0) next
    circos.text(
      x = n_in_sec / 2, y = 0.5,
      labels = if (sec == "Singletons") "Singletons" else paste0("C", sub("^Cluster ", "", sec)),
      sector.index = sec, track.index = clust_ti,
      facing = "bending.inside", niceFacing = TRUE,
      cex = 0.25, col = "white", adj = c(0.5, 0.5)
    )
  }

  circos.text(
    x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = clust_ti)[2],
    y = 0.5, labels = "Instrument Cluster ",
    sector.index = last_sector, track.index = clust_ti,
    facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
  )

  # ── Forward group label track ─────────────────────────────────────────────
  circos.track(
    ylim = c(0, 1), track.height = 0.018, track.margin = c(0.001, 0.001),
    bg.col = "white", bg.border = NA,
    panel.fun = function(x, y) {}
  )
  fwd_label_ti <- get.current.track.index()
  op <- par(xpd = NA)
  circos.text(
    x = 90 + 27, y = 0.5,
    labels = paste0("Metabolite → ", this_outcome, " MR"),
    sector.index = first_sec, track.index = fwd_label_ti,
    facing = "bending.inside", niceFacing = TRUE,
    cex = 0.4, font = 2, col = "grey20", adj = c(0.5, 0.5)
  )
  par(op)

  # ── Forward nSNP track ────────────────────────────────────────────────────
  circos.track(
    ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
    bg.border = NA, bg.col = NA,
    panel.fun = function(x, y) {
      sec  <- get.cell.meta.data("sector.index")
      mets <- sector_mets[[sec]]
      vals <- nsnp_fwd_i[match(mets, metab_order)]
      for (k in seq_along(mets))
        if (!is.na(vals[k]))
          circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                      facing = "clockwise", niceFacing = TRUE,
                      cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
    }
  )
  nsnp_fwd_ti <- get.current.track.index()
  circos.text(
    x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = nsnp_fwd_ti)[2],
    y = 0.5, labels = "Fwd nSNP ",
    sector.index = last_sector, track.index = nsnp_fwd_ti,
    facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
  )

  # ── Forward Fstat track ───────────────────────────────────────────────────
  circos.track(
    ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
    bg.border = NA, bg.col = NA,
    panel.fun = function(x, y) {
      sec  <- get.cell.meta.data("sector.index")
      mets <- sector_mets[[sec]]
      vals <- fstat_fwd_i[match(mets, metab_order)]
      for (k in seq_along(mets))
        if (!is.na(vals[k]))
          circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                      facing = "clockwise", niceFacing = TRUE,
                      cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
    }
  )
  fstat_fwd_ti <- get.current.track.index()
  circos.text(
    x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = fstat_fwd_ti)[2],
    y = 0.5, labels = "Fwd Fstat ",
    sector.index = last_sector, track.index = fstat_fwd_ti,
    facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
  )

  # ── Forward Steiger-excluded nSNP track ───────────────────────────────────
  steiger_excl_fwd_i <- nsnp_fwd_i - steiger_nsnp_fwd_i
  circos.track(
    ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
    bg.border = NA, bg.col = NA,
    panel.fun = function(x, y) {
      sec  <- get.cell.meta.data("sector.index")
      mets <- sector_mets[[sec]]
      vals <- steiger_excl_fwd_i[match(mets, metab_order)]
      for (k in seq_along(mets))
        if (!is.na(vals[k]))
          circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                      facing = "clockwise", niceFacing = TRUE,
                      cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
    }
  )
  steiger_nsnp_fwd_ti <- get.current.track.index()
  circos.text(
    x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = steiger_nsnp_fwd_ti)[2],
    y = 0.5, labels = "Steiger excl. ",
    sector.index = last_sector, track.index = steiger_nsnp_fwd_ti,
    facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
  )

  # ── Forward MR heatmap ────────────────────────────────────────────────────
  circos.heatmap(
    circ_mat_fwd, split = pathway_split, col = col_disc_i, na.col = "white",
    cluster = FALSE, cell.border = "gray95", cell.lwd = 0.5,
    track.height = 0.13, show.sector.labels = FALSE, rownames.side = "none"
  )
  fwd_ti <- get.current.track.index()

  # Egger * overlay
  sig_map_i <- egger_intercept_sig[[this_outcome]]
  for (sec in names(sector_mets)) {
    mets <- sector_mets[[sec]]
    for (k in seq_along(mets))
      if (isTRUE(sig_map_i[mets[k]]))
        circos.text(k - 0.5, egger_y_per_out, labels = "*",
                    sector.index = sec, track.index = fwd_ti,
                    facing = "inside", niceFacing = TRUE,
                    cex = 0.4, col = "grey40", adj = c(0.5, 0.5))
  }

  # Method labels
  x_right_fwd <- get.cell.meta.data("xlim", sector.index = last_sector, track.index = fwd_ti)[2]
  for (j in seq_len(n_meth_per_out))
    circos.text(x_right_fwd, j - 0.5, labels = paste0(meth_abbr_per_out[j], " "),
                sector.index = last_sector, track.index = fwd_ti,
                facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5))

  # ── Reverse MR heatmap ────────────────────────────────────────────────────
  if (!is.null(rev_mat)) {
    # ── Reverse group label track ───────────────────────────────────────────
    circos.track(
      ylim = c(0, 1), track.height = 0.018, track.margin = c(0.001, 0.001),
      bg.col = "white", bg.border = NA,
      panel.fun = function(x, y) {}
    )
    rev_label_ti <- get.current.track.index()
    op <- par(xpd = NA)
    circos.text(
      x = 90 + 28, y = 0.5,
      labels = paste0(this_outcome, " → Metabolite MR"),
      sector.index = first_sec, track.index = rev_label_ti,
      facing = "bending.inside", niceFacing = TRUE,
      cex = 0.4, font = 2, col = "grey20", adj = c(0.5, 0.5)
    )
    par(op)

    # ── Reverse nSNP track ────────────────────────────────────────────────
    circos.track(
      ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
      bg.border = NA, bg.col = NA,
      panel.fun = function(x, y) {
        sec  <- get.cell.meta.data("sector.index")
        mets <- sector_mets[[sec]]
        vals <- nsnp_rev_i[match(mets, metab_order)]
        for (k in seq_along(mets))
          if (!is.na(vals[k]))
            circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                        facing = "clockwise", niceFacing = TRUE,
                        cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
      }
    )
    nsnp_rev_ti <- get.current.track.index()
    circos.text(
      x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = nsnp_rev_ti)[2],
      y = 0.5, labels = "Rev nSNP ",
      sector.index = last_sector, track.index = nsnp_rev_ti,
      facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
    )

    # ── Reverse Fstat track ───────────────────────────────────────────────
    circos.track(
      ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
      bg.border = NA, bg.col = NA,
      panel.fun = function(x, y) {
        sec  <- get.cell.meta.data("sector.index")
        mets <- sector_mets[[sec]]
        vals <- fstat_rev_i[match(mets, metab_order)]
        for (k in seq_along(mets))
          if (!is.na(vals[k]))
            circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                        facing = "clockwise", niceFacing = TRUE,
                        cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
      }
    )
    fstat_rev_ti <- get.current.track.index()
    circos.text(
      x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = fstat_rev_ti)[2],
      y = 0.5, labels = "Rev Fstat ",
      sector.index = last_sector, track.index = fstat_rev_ti,
      facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
    )

    # ── Reverse Steiger-excluded nSNP track ───────────────────────────────
    steiger_excl_rev_i <- nsnp_rev_i - steiger_nsnp_rev_i
    circos.track(
      ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
      bg.border = NA, bg.col = NA,
      panel.fun = function(x, y) {
        sec  <- get.cell.meta.data("sector.index")
        mets <- sector_mets[[sec]]
        vals <- steiger_excl_rev_i[match(mets, metab_order)]
        for (k in seq_along(mets))
          if (!is.na(vals[k]))
            circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                        facing = "clockwise", niceFacing = TRUE,
                        cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
      }
    )
    steiger_nsnp_rev_ti <- get.current.track.index()
    circos.text(
      x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = steiger_nsnp_rev_ti)[2],
      y = 0.5, labels = "Steiger excl. ",
      sector.index = last_sector, track.index = steiger_nsnp_rev_ti,
      facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
    )

    circos.heatmap(
      rev_mat, split = pathway_split, col = col_disc_i, na.col = "white",
      cluster = FALSE, cell.border = "gray95", cell.lwd = 0.5,
      track.height = 0.13, show.sector.labels = FALSE, rownames.side = "none"
    )
    rev_ti <- get.current.track.index()

    x_right_rev <- get.cell.meta.data("xlim", sector.index = last_sector, track.index = rev_ti)[2]
    for (j in seq_len(n_meth_per_out))
      circos.text(x_right_rev, j - 0.5, labels = paste0(meth_abbr_per_out[j], " "),
                  sector.index = last_sector, track.index = rev_ti,
                  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5))

  }

  legend("topright", inset = c(0.1, 0.15),
         legend = c("FDR p<0.05 (risk)", "p<0.05 (risk)", "p≥05 (risk)",
                    "p≥05 (protective)", "p<0.05 (protective)", "FDR p<0.05 (protective)"),
         fill   = c(c_red_fdr, c_red_nom, c_red_ns, c_blue_ns, c_blue_nom, c_blue_fdr),
         title  = "MR result (FDR within outcome × method)", bty = "n", cex = 0.55)

  circos.clear()
  dev.off()
  cat("Saved:", out_circos_file, "\n")
}
