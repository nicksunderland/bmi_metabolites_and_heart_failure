#### DESCRIPTION #######################################
# Circos plot of BMI → metabolite MR + observational results
########################################################

#### SET INPUT #########################################
results_file     <- snakemake@input[["results_file"]]
name_map_file    <- snakemake@input[["name_map_file"]]
replicating_file <- snakemake@input[["replicating_file"]]
circos_file      <- snakemake@output[["circos"]]
scatter_file     <- snakemake@output[["scatter"]]
log_file         <- snakemake@log[["log"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
results_file     <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/bmi_on_metabolites.tsv")
name_map_file    <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "scripts/gwas_metab_name_map.xlsx")
replicating_file <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/replication/study_replication.tsv")
circos_file      <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/figures/bmi_mr/bmi_metab_mr_circos.png")
scatter_file     <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/figures/bmi_mr/bmi_metab_mr_scatter.png")
log_file         <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/logs/bmi_metab_mr_discordance.log")
# ─────────────────────────────────────────────────────────────────────────────
}

library(data.table)
library(readxl)
library(circlize)
library(ggplot2)

dt  <- fread(results_file, sep = "\t")
dt_steiger <- dt[steiger_filtering == TRUE & method == "Inverse variance weighted"]
dt_steiger[, method := paste(method, "Steiger filtered")]
dt  <- rbind(
  dt[steiger_filtering == FALSE],
  dt_steiger
)[order(exp_id, out_id, method)]
rep <- fread(replicating_file, na.strings = "-")
map <- read_xlsx(name_map_file, sheet = 1) |> as.data.table()

# Convert observational columns stored as strings
obs_cols <- intersect(
  c("estimate_bbs", "estimate_direct", "std.error_bbs", "std.error_direct",
    "p.value_bbs", "p.value_fdr_bbs",
    "p.value_interaction_direct", "p.value_interaction_fdr_direct"),
  names(rep)
)
if (length(obs_cols)) rep[, (obs_cols) := lapply(.SD, as.numeric), .SDcols = obs_cols]

# Join pathway info
dt[map, c("label", "pathway_group", "platform") := .(i.label, i.pathway_group, i.platform),
   on = c("out_id" = "gwas_id")]
dt[, pathway_broad := fcase(
  pathway_group == "Lipid metabolism",      "Lipid metabolism",
  pathway_group == "Lipoproteins",          "Lipoproteins",
  pathway_group == "Amino acid metabolism", "Amino acids",
  default = "Other"
)]
dt[, pathway_broad := factor(pathway_broad,
                             levels = c("Lipid metabolism", "Lipoproteins", "Amino acids", "Other"))]
dt[, label := fcoalesce(label, out_id)]

# Save pre-disambiguation label for joining rep (which uses original label + platform)
dt[, raw_label := label]

# Filter to Trial consistent
dt[rep, replicates := i.replicates, on = c("raw_label" = "label", "platform")]
dt <- dt[replicates == "Trial consistent (interaction)"]

# Disambiguate labels that appear on more than one platform
platform_suffix <- c("nightingale" = "NMR", "metabolon" = "MS")
dup_labels <- dt[, .(n = uniqueN(out_id)), by = label][n > 1, label]
dt[label %in% dup_labels, label := paste0(label, " (", platform_suffix[platform], ")")]

# Subset to circos MR methods
circ_method_levels <- c(
  "Inverse variance weighted",
  "Inverse variance weighted Steiger filtered",
  # "MR-PRESSO-corrected", # doesnt work well for very large instruments as no outliers found, other sensitivity analyses reasonable here
  "Weighted median",
  "Weighted mode",
  "MR Egger"
)
dt <- dt[method %in% circ_method_levels]
dt[, method := factor(method, levels = circ_method_levels)]

# ── Metabolite ordering: by pathway then BBS effect size (largest positive first) ──
bbs_ord <- unique(dt[method == "Inverse variance weighted",
                     .(label = as.character(label), raw_label, platform, pathway_broad)])
bbs_ord[rep, bbs_b := i.estimate_bbs, on = c("raw_label" = "label", "platform")]
bbs_ord <- bbs_ord[, .(bbs_b = mean(bbs_b, na.rm = TRUE), pathway_broad = pathway_broad[1L]),
                   by = label]
get_lvl <- function(grp)
  bbs_ord[pathway_broad == grp][order(-bbs_b, na.last = TRUE), as.character(label)]
circ_metab_order <- c(
  get_lvl("Lipid metabolism"),
  get_lvl("Lipoproteins"),
  get_lvl("Amino acids"),
  get_lvl("Other")
)
dt[, label := factor(label, levels = circ_metab_order)]

pw_lookup     <- unique(dt[, .(label = as.character(label), pathway_broad)])
pathway_split <- factor(
  pw_lookup[match(circ_metab_order, label), as.character(pathway_broad)],
  levels = c("Lipid metabolism", "Lipoproteins", "Amino acids", "Other")
)

n_meth_circ <- length(circ_method_levels)
last_sector  <- levels(pathway_split)[length(levels(pathway_split))]
first_sec    <- levels(pathway_split)[1]
sector_mets  <- split(circ_metab_order, pathway_split)

# ── Colour scheme for MR heatmap (direction × significance tiers) ────────────
col_blue <- "#2166AC"; col_red <- "#B2182B"
c_blue_fdr <- adjustcolor(col_blue, alpha.f = 1.00)
c_blue_nom <- adjustcolor(col_blue, alpha.f = 0.50)
c_blue_ns  <- adjustcolor(col_blue, alpha.f = 0.15)
c_red_fdr  <- adjustcolor(col_red,  alpha.f = 1.00)
c_red_nom  <- adjustcolor(col_red,  alpha.f = 0.50)
c_red_ns   <- adjustcolor(col_red,  alpha.f = 0.15)

sig_colour <- function(b, p, p_fdr) {
  fcase(
    is.na(p) | is.na(b),  "white",
    b < 0 & p_fdr < 0.05, c_blue_fdr,
    b < 0 & p < 0.05,     c_blue_nom,
    b < 0,                 c_blue_ns,
    b > 0 & p_fdr < 0.05, c_red_fdr,
    b > 0 & p < 0.05,     c_red_nom,
    b > 0,                 c_red_ns,
    default               = "white"
  )
}

# ── MR colour matrix ──────────────────────────────────────────────────────────
dt[, p_fdr := p.adjust(p, method = "fdr"), by = method]
dt[, col   := sig_colour(b, p, p_fdr)]

# ── Direction discordance log ─────────────────────────────────────────────────
{
  # IVW rows with valid estimates; join trial direction (BBS, falling back to DiRECT)
  dt_ivw <- dt[method == "Inverse variance weighted" & !is.na(b) & !is.na(p)]
  dt_ivw[rep,
         trial_b := fcoalesce(as.numeric(i.estimate_bbs), as.numeric(i.estimate_direct)),
         on = c("raw_label" = "label", "platform")]

  # a) metabolites with MR data
  n_mr <- uniqueN(dt_ivw$label)

  # b/c) discordant: IVW FDR p<0.05 AND sign flip vs trial
  disc_labels <- dt_ivw[
    p_fdr < 0.05 & !is.na(trial_b) & sign(b) != sign(trial_b),
    unique(as.character(label))
  ]
  n_disc   <- length(disc_labels)
  pct_disc <- round(100 * n_disc / n_mr, 1)

  # d/e) discordant Steiger: Steiger-filtered IVW FDR p<0.05 AND flipped vs trial
  dt_st <- dt[method == "Inverse variance weighted Steiger filtered" & !is.na(b)]
  dt_st[rep,
        trial_b := fcoalesce(as.numeric(i.estimate_bbs), as.numeric(i.estimate_direct)),
        on = c("raw_label" = "label", "platform")]
  disc_st_labels <- dt_st[
    p_fdr < 0.05 & !is.na(trial_b) & sign(b) != sign(trial_b),
    unique(as.character(label))
  ]
  n_disc_st   <- length(disc_st_labels)
  pct_disc_st <- round(100 * n_disc_st / n_mr, 1)

  # f/g) strongly discordant: discordant IVW AND discordant Steiger
  strong_labels <- intersect(disc_labels, disc_st_labels)
  n_strong      <- length(strong_labels)
  pct_strong    <- round(100 * n_strong / n_mr, 1)

  log_lines <- c(
    "=== BMI -> Metabolite MR: Direction Discordance vs Trial Estimates ===",
    "",
    sprintf("a) Metabolites with MR data available:                        %d", n_mr),
    sprintf("b) Discordant IVW (IVW FDR p<0.05 + sign flip):              %d", n_disc),
    sprintf("c) Discordant IVW %%:                                         %.1f%%", pct_disc),
    sprintf("d) Discordant Steiger (Steiger FDR p<0.05 + sign flip):       %d", n_disc_st),
    sprintf("e) Discordant Steiger %%:                                      %.1f%%", pct_disc_st),
    sprintf("f) Strongly discordant (discordant on both IVW + Steiger):    %d", n_strong),
    sprintf("g) Strongly discordant %%:                                     %.1f%%", pct_strong),
    "",
    "--- Definitions ---",
    "Discordant IVW:     IVW FDR p<0.05 AND sign(IVW beta) != sign(trial beta)",
    "Discordant Steiger: Steiger-filtered IVW FDR p<0.05 AND sign(Steiger beta) != sign(trial beta)",
    "Strongly discordant: meets both Discordant IVW AND Discordant Steiger criteria",
    "Trial direction:    sign(estimate_bbs), falling back to sign(estimate_direct) if missing"
  )

  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(log_lines, log_file)
}

wide_mr <- dcast(dt, label ~ method, value.var = "col", fun.aggregate = function(x) x[1L])
mr_mat  <- as.matrix(wide_mr[match(circ_metab_order, as.character(wide_mr$label)),
                              circ_method_levels, with = FALSE])
mr_mat[is.na(mr_mat)] <- "white"
rownames(mr_mat) <- circ_metab_order
col_disc <- setNames(unique(as.vector(mr_mat)), unique(as.vector(mr_mat)))

# ── Egger intercept flag (p < 0.05 → overlay "*" on Egger column) ────────────
egger_p_dt <- dt[as.character(method) == "MR Egger",
                 .(label = as.character(label), egger_intercept_p)]
egger_sig  <- setNames(
  !is.na(egger_p_dt$egger_intercept_p) & egger_p_dt$egger_intercept_p < 0.05,
  egger_p_dt$label
)
# column j is at y = n_meth_circ - j + 0.5 in the heatmap coordinate system
egger_y    <- n_meth_circ - match("MR Egger", circ_method_levels) + 0.5

# ── Per-metabolite instrument info (from IVW row) ────────────────────────────
dt[, nsnp_steiger_filtered := nsnp_steiger_filtered[method=="Inverse variance weighted Steiger filtered"], by = "out_id"]
ivw_info <- dt[method == "Inverse variance weighted",
               .(nsnp         = nsnp[1L],
                 fstat        = if ("fstat"               %in% names(dt)) fstat[1L]               else NA_real_,
                 steiger_snps = nsnp_steiger_filtered[1L]),
               by = label]
nsnp_vec    <- ivw_info[match(circ_metab_order, as.character(label)), nsnp]
fstat_vec   <- ivw_info[match(circ_metab_order, as.character(label)), fstat]
steiger_vec <- ivw_info[match(circ_metab_order, as.character(label)), steiger_snps]

# ── Observational estimates: join rep via original label + platform ───────────
metab_meta <- unique(dt[, .(label = as.character(label), raw_label, platform)])[
  match(circ_metab_order, label)]
rep_ord <- rep[metab_meta, .(bbs_b = estimate_bbs, dir_b = estimate_direct),
               on = c("label" = "raw_label", "platform")]
bbs_vec <- rep_ord$bbs_b
dir_vec <- rep_ord$dir_b

bbs_col     <- ifelse(is.na(bbs_vec), "white", ifelse(bbs_vec > 0, c_red_fdr, c_blue_fdr))
dir_col     <- ifelse(is.na(dir_vec), "white", ifelse(dir_vec > 0, c_red_fdr, c_blue_fdr))
bbs_mat     <- matrix(bbs_col, ncol = 1, dimnames = list(circ_metab_order, ""))
dir_mat     <- matrix(dir_col, ncol = 1, dimnames = list(circ_metab_order, ""))
obs_col_map <- setNames(unique(c(bbs_col, dir_col)), unique(c(bbs_col, dir_col)))

# ── Draw circos ───────────────────────────────────────────────────────────────
meth_abbr    <- rev(c("IVW ", "IVW Steiger ", "W.Med ", "W.Mode ", "Egger "))
gap          <- 18
label_sec    <- levels(pathway_split)[2]   # sector over which ring labels are placed; change index to reposition
label_x_frac <- 0.74                        # 0 = sector start, 1 = sector end; fine-tune within that sector
label_y_mr   <- 5.5                        # y for MR ring label; > n_meth_circ (5) = outer margin, tune downward
label_y_obs  <- 1.4                        # y for DiRECT/BBS ring labels; > 1 = outer margin, tune upward
dir.create(dirname(circos_file), recursive = TRUE, showWarnings = FALSE)
png(circos_file, width = 10, height = 10, units = "in", res = 600, bg = "white")
circos.clear()
circos.par(
  gap.after               = setNames(c(3, 3, 3, gap), levels(pathway_split)),
  track.margin            = c(0.005, 0.018),
  cell.padding            = c(0, 0, 0, 0),
  start.degree            = 270 - gap,
  points.overflow.warning = FALSE,
  canvas.xlim             = c(-1.1, 1.1),
  canvas.ylim             = c(-1.1, 1.1)
)

# Track 1 (outermost): invisible thin heatmap — initialises circle, labels appear outside
label_mat <- matrix(rep("white", length(circ_metab_order)), ncol = 1,
                    dimnames = list(circ_metab_order, ""))
circos.heatmap(
  label_mat,
  split              = pathway_split,
  col                = c("white" = "white"),
  cluster            = FALSE,
  track.height       = 0.001,
  track.margin       = c(0.002, 0.002),
  show.sector.labels = FALSE,
  rownames.side      = "outside",
  rownames.cex       = 0.28,
  rownames.font      = 1,
  cell.border        = NA,
  na.col             = "white"
)

# Track 2: nSNPs
circos.track(
  ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
  bg.border = NA, bg.col = NA,
  panel.fun = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    mets <- sector_mets[[sec]]
    vals <- nsnp_vec[match(mets, circ_metab_order)]
    for (k in seq_along(mets))
      if (!is.na(vals[k]))
        circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
  }
)
nsnp_tidx <- get.current.track.index()
circos.text(
  x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = nsnp_tidx)[2],
  y = 0.5, labels = "nSNPs ",
  sector.index = last_sector, track.index = nsnp_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)

# Track 3: Fstat
circos.track(
  ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
  bg.border = NA, bg.col = NA,
  panel.fun = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    mets <- sector_mets[[sec]]
    vals <- fstat_vec[match(mets, circ_metab_order)]
    for (k in seq_along(mets))
      if (!is.na(vals[k]))
        circos.text(k - 0.5, 0.5, labels = round(vals[k], 0),
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
  }
)
fstat_tidx <- get.current.track.index()
circos.text(
  x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = fstat_tidx)[2],
  y = 0.5, labels = "Fstat ",
  sector.index = last_sector, track.index = fstat_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)

# Track 4: Steiger filtered SNPs
circos.track(
  ylim = c(0, 1), track.height = 0.012, track.margin = c(0.002, 0.002),
  bg.border = NA, bg.col = NA,
  panel.fun = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    mets <- sector_mets[[sec]]
    vals <- steiger_vec[match(mets, circ_metab_order)]
    for (k in seq_along(mets)) {
      if (!is.na(vals[k])) {
        circos.text(k - 0.5, 0.5, labels = vals[k],
                    facing = "clockwise", niceFacing = TRUE,
                    cex = 0.2, col = "grey20", adj = c(0.5, 0.5))
      }
    }
  }
)
stei_tidx <- get.current.track.index()
circos.text(
  x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = stei_tidx)[2],
  y = 0.5, labels = "Steiger excluded nSNPs ",
  sector.index = last_sector, track.index = stei_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)

# Track 5: BMI → Metabolite MR heatmap (multi-method columns)
circos.heatmap(
  mr_mat,
  split              = pathway_split,
  col                = col_disc,
  na.col             = "white",
  cluster            = FALSE,
  cell.border        = "gray95",
  cell.lwd           = 0.5,
  track.height       = 0.11,
  show.sector.labels = FALSE,
  rownames.side      = "none"
)
mr_tidx    <- get.current.track.index()
x_right_mr <- get.cell.meta.data("xlim", sector.index = last_sector, track.index = mr_tidx)[2]
for (j in seq_len(n_meth_circ))
  circos.text(x = x_right_mr, y = j - 0.5, labels = meth_abbr[j],
              sector.index = last_sector, track.index = mr_tidx,
              facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5))
for (sec in names(sector_mets)) {
  mets <- sector_mets[[sec]]
  for (k in seq_along(mets))
    if (isTRUE(egger_sig[mets[k]]))
      circos.text(k - 0.5, egger_y, labels = "*",
                  sector.index = sec, track.index = mr_tidx,
                  facing = "inside", niceFacing = TRUE,
                  cex = 0.4, col = "grey40", adj = c(0.5, 0.5))
}
{
  top_r   <- get.cell.meta.data("cell.top.radius", sector.index = label_sec, track.index = mr_tidx)
  draw.sector(start.degree = 0, end.degree = 360,
              rou1 = top_r + 0.017, rou2 = top_r + 0.001, col = "white", border = NA)
  xlim_mr <- get.cell.meta.data("xlim", sector.index = label_sec, track.index = mr_tidx)
  op <- par(xpd = NA)
  circos.text(
    x            = xlim_mr[1] + label_x_frac * diff(xlim_mr),
    y            = label_y_mr,
    labels       = "Life-time BMI exposure effect (MR)",
    sector.index = label_sec,
    track.index  = mr_tidx,
    facing       = "bending.inside",
    niceFacing   = TRUE,
    cex          = 0.4,
    font         = 2,
    col          = "grey20",
    adj          = c(0.5, 0.5)
  )
  par(op)
}

# Track 6: DiRECT observational estimate (gradient: blue=negative, red=positive)
circos.heatmap(
  dir_mat,
  split              = pathway_split,
  col                = obs_col_map,
  na.col             = "white",
  cluster            = FALSE,
  cell.border        = NA,
  track.height       = 0.03,
  show.sector.labels = FALSE,
  rownames.side      = "none"
)
dir_tidx <- get.current.track.index()
circos.text(
  x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = dir_tidx)[2],
  y = 0.5, labels = "DiRECT ",
  sector.index = last_sector, track.index = dir_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)
{
  dir_top_r <- get.cell.meta.data("cell.top.radius", sector.index = label_sec, track.index = dir_tidx)
  draw.sector(start.degree = 0, end.degree = 360,
              rou1 = dir_top_r + 0.008, rou2 = dir_top_r + 0.001, col = "white", border = NA)
  xlim_dir  <- get.cell.meta.data("xlim", sector.index = label_sec, track.index = dir_tidx)
  op <- par(xpd = NA)
  circos.text(
    x            = xlim_dir[1] + label_x_frac * diff(xlim_dir),
    y            = label_y_obs,
    labels       = "Structured dietary programme effect",
    sector.index = label_sec,
    track.index  = dir_tidx,
    facing       = "bending.inside",
    niceFacing   = TRUE,
    cex          = 0.4,
    font         = 2,
    col          = "grey20",
    adj          = c(0.5, 0.5)
  )
  par(op)
}

# Track 7: BBS observational estimate
circos.heatmap(
  bbs_mat,
  split              = pathway_split,
  col                = obs_col_map,
  na.col             = "white",
  cluster            = FALSE,
  cell.border        = NA,
  track.height       = 0.03,
  show.sector.labels = FALSE,
  rownames.side      = "none"
)
bbs_tidx <- get.current.track.index()
circos.text(
  x = get.cell.meta.data("xlim", sector.index = last_sector, track.index = bbs_tidx)[2],
  y = 0.5, labels = "BBS ",
  sector.index = last_sector, track.index = bbs_tidx,
  facing = "inside", niceFacing = TRUE, cex = 0.3, adj = c(0, 0.5)
)
{
  bbs_top_r <- get.cell.meta.data("cell.top.radius", sector.index = label_sec, track.index = bbs_tidx)
  draw.sector(start.degree = 0, end.degree = 360,
              rou1 = bbs_top_r + 0.008, rou2 = bbs_top_r + 0.001, col = "white", border = NA)
  xlim_bbs  <- get.cell.meta.data("xlim", sector.index = label_sec, track.index = bbs_tidx)
  op <- par(xpd = NA)
  circos.text(
    x            = xlim_bbs[1] + label_x_frac * diff(xlim_bbs),
    y            = label_y_obs,
    labels       = "Bariatric surgery effect",
    sector.index = label_sec,
    track.index  = bbs_tidx,
    facing       = "bending.inside",
    niceFacing   = TRUE,
    cex          = 0.4,
    font         = 2,
    col          = "grey20",
    adj          = c(0.5, 0.5)
  )
  par(op)
}

# Track 8 (innermost): pathway group coloured band — large white centre left free
pw_cols <- c(
  "Lipid metabolism" = "#E8A838",
  "Lipoproteins"     = "#6BAED6",
  "Amino acids"      = "#74C476",
  "Other"            = "#BCBCBC"
)
circos.track(
  ylim = c(0, 1), track.height = 0.03, track.margin = c(0.005, 0.005),
  bg.border = NA,
  panel.fun = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = pw_cols[sec], border = NA)
    circos.text(mean(xlim), mean(ylim), labels = sec,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.45, col = "white", font = 2)
  }
)

legend("topright", inset = 0.15,
       legend = c("FDR p<0.05 (↑ effect)", "p<0.05 (↑ effect)", "n.s. (↑ effect)",
                  "n.s. (↓ effect)", "p<0.05 (↓ effect)", "FDR p<0.05 (↓ effect)"),
       fill   = c(c_red_fdr, c_red_nom, c_red_ns, c_blue_ns, c_blue_nom, c_blue_fdr),
       title  = "Exposure effect on metabolite\n(FDR within intervention / MR method)", bty = "n", cex = 0.5)

circos.clear()
dev.off()
cat("Saved:", circos_file, "\n")



# ── Scatter: trial estimate (x) vs BMI → metabolite MR estimate (y) ──────────
# Negative slope expected if BMI mediates trial effects on metabolites
ivw_mr <- dt[method == "Inverse variance weighted",
             .(label = as.character(label), mr_b = b, mr_p_fdr = p_fdr)]

scatter_bbs <- rep[metab_meta,
                   .(label = i.label, trial_b = estimate_bbs, trial_fdr = p.value_fdr_bbs),
                   on = c("label" = "raw_label", "platform")]
scatter_bbs[, trial := "Bariatric surgery (BBS)"]

scatter_dir <- rep[metab_meta,
                   .(label = i.label, trial_b = estimate_direct,
                     trial_fdr = p.value_interaction_fdr_direct),
                   on = c("label" = "raw_label", "platform")]
scatter_dir[, trial := "Dietary programme (DiRECT)"]

scatter_dt <- rbindlist(list(scatter_bbs, scatter_dir), use.names = TRUE)
scatter_dt[ivw_mr, c("mr_b", "mr_p_fdr") := .(i.mr_b, i.mr_p_fdr), on = "label"]
scatter_dt <- scatter_dt[!is.na(trial_b) & !is.na(mr_b)]
trial_cols  <- c("Bariatric surgery (BBS)" = "#B2182B", "Dietary programme (DiRECT)" = "#2166AC")
trial_shape <- c("Bariatric surgery (BBS)" = 16L,       "Dietary programme (DiRECT)" = 17L)

p_scatter <- ggplot(scatter_dt,
                    aes(x = trial_b, y = mr_b, colour = trial, shape = trial)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.12,
              formula = y ~ x) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = trial_cols, name = NULL) +
  scale_shape_manual(values  = trial_shape, name = NULL) +
  labs(
    x = "Weight loss trial effect on metabolite\n(SD change with intervention)",
    y = "Life-time BMI effect on metabolite - IVW MR\n(SD change per 1-SD higher BMI)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background    = element_blank(),
    legend.key.size      = unit(0.4, "cm"),
    legend.text          = element_text(size = 8),
    axis.line            = element_line(colour = "grey50", linewidth = 0.4),
    axis.ticks           = element_line(colour = "grey50", linewidth = 0.3),
    axis.text            = element_text(colour = "grey30", size = 8),
    axis.title           = element_text(colour = "grey20", size = 12),
    panel.grid           = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)
  )

ggsave(scatter_file, p_scatter, width = 5, height = 5, dpi = 300, bg = "transparent")
cat("Saved:", scatter_file, "\n")
