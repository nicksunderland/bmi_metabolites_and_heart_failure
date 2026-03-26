#### DESCRIPTION #######################################
# Per-pathway cis-MR forest plots

#### SET INPUT #########################################
mr_files    <- snakemake@input[["mr_files"]]
gw_mr_file  <- snakemake@input[["gw_mr"]]
pathway     <- snakemake@params[["pathway"]]
gwas_ids    <- snakemake@params[["gwas_ids"]]
forest_plot <- snakemake@output[["forest_plot"]]
########################################################

message("[DEBUG] script start")
message("[DEBUG] pathway:     ", pathway)
message("[DEBUG] forest_plot: ", forest_plot)
message("[DEBUG] mr_files length: ", length(mr_files))
for (f in mr_files) message("[DEBUG]   file: ", f, "  exists=", file.exists(f))
message("[DEBUG] gw_mr_file: ", gw_mr_file, "  exists=", file.exists(gw_mr_file))
message("[DEBUG] gwas_ids: ", paste(names(gwas_ids), unlist(gwas_ids), sep = "=", collapse = ", "))

dir.create(dirname(forest_plot), recursive = TRUE, showWarnings = FALSE)

# testing
if (FALSE) {
  library(data.table)
  pathway    <- "branched_chain_amino_acids"
  repo_dir   <- Sys.getenv("HF_METABOLITE_REPO2")
  mr_files   <- list.files(
    file.path(repo_dir, "output", "tmp_objects", "pathway_cis_mr", pathway),
    pattern = "\\.tsv\\.gz$", full.names = TRUE
  )
  gw_mr_file <- file.path(repo_dir, "output", "tables", "mr_results",
                           "metabolite_outcome_mr_results.tsv.gz")
  gwas_ids   <- list(valine    = "GCST90302122",
                     leucine   = "GCST90301994",
                     isoleucine = "GCST90301987",
                     kiv       = "GCST90199677",
                     kic       = "GCST90199655",
                     kmv       = "GCST90199631")
  forest_plot <- file.path(repo_dir, "output", "figures", "cis_mr",
                           paste0(pathway, "_mr_forest.png"))
}

message("[DEBUG] loading libraries...")
library(data.table)
library(ggplot2)
message("[DEBUG] libraries loaded")

fmt_name <- function(x) {
  x <- trimws(gsub("_", " ", x))
  ifelse(nchar(x) <= 4, toupper(x), tools::toTitleCase(x))
}

outcome_order  <- c("heart_failure", "hfref", "hfpef")
outcome_labels <- c(heart_failure = "All HF", hfref = "HFrEF", hfpef = "HFpEF")

method_order   <- c("mr_ivw", "mr_egger", "mr_weighted_median", "mr_weighted_mode")
method_labels  <- c(
  mr_ivw             = "IVW",
  mr_egger           = "MR-Egger",
  mr_weighted_median = "Weighted median",
  mr_weighted_mode   = "Weighted mode"
)
method_colours <- c(
  "IVW"             = "#1b4f8a",
  "MR-Egger"        = "#c0392b",
  "Weighted median" = "#229954",
  "Weighted mode"   = "#8e44ad"
)
method_shapes <- c(
  "IVW"             = 16,
  "MR-Egger"        = 17,
  "Weighted median" = 15,
  "Weighted mode"   = 18
)
method_sizes <- c(
  "IVW"             = 2.8,
  "MR-Egger"        = 2.1,
  "Weighted median" = 2.1,
  "Weighted mode"   = 2.1
)

# read cis-MR
message("[DEBUG] reading ", length(mr_files), " cis-MR result files")
dat <- rbindlist(lapply(mr_files, fread), fill = TRUE)
message("[DEBUG] cis read: ", nrow(dat), " rows; cols: ", paste(names(dat), collapse = ", "))

dat_cis     <- dat[is.na(error)]
dat_cis_err <- dat[!is.na(error)]
if (nrow(dat_cis_err) > 0) {
  message("[warn] ", nrow(dat_cis_err), " cis error rows excluded:")
  message(paste(dat_cis_err[, paste(gwas_name, gene_name, outcome, error, sep = " | ")],
                collapse = "\n"))
}
message("[DEBUG] dat_cis: ", nrow(dat_cis), " ok rows; ", nrow(dat_cis_err), " error rows")

# read genome-wide MR
message("[DEBUG] reading genome-wide MR results: ", gw_mr_file)
gw_all <- fread(gw_mr_file)
message("[DEBUG] gw_all: ", nrow(gw_all), " rows",
        if (ncol(gw_all) > 0) paste0("; cols: ", paste(names(gw_all), collapse = ", ")) else "")

# lookup
id_map <- setNames(names(gwas_ids), unlist(gwas_ids))
id_map <- c(id_map, setNames(names(gwas_ids), paste0(unlist(gwas_ids), ".fst")))
id_map <- c(id_map, setNames(names(gwas_ids), names(gwas_ids)))
match_col <- if ("metabolite_id" %in% names(gw_all)) "metabolite_id" else "exposure"
message("[DEBUG] matching gw_all on column: ", match_col)

gw_sub <- gw_all[get(match_col) %in% names(id_map) & (is.na(error)|error=="")]
message("[DEBUG] gw_sub after ID filter: ", nrow(gw_sub), " rows")

if (nrow(gw_sub) > 0) {
  gw_sub[, gwas_name := id_map[get(match_col)]]
  gw_sub[, gene_name := "Genome-wide"]
  gw_sub[, pathway   := pathway]
} else {
  message("[warn] no genome-wide rows matched for pathway: ", pathway,
          " (", match_col, " values in file: ",
          paste(head(unique(gw_all[[match_col]]), 5), collapse = ", "), "...)")
}

# combine
common_cols <- c("exposure", "gwas_name", "gene_name", "pathway",
                 "outcome", "method", "n_snp", "fstat",
                 "b", "b_se", "p", "error")

ensure_cols <- function(dt, cols) {
  for (col in setdiff(cols, names(dt))) dt[, (col) := NA]
  dt[, ..cols]
}

dat_cis <- ensure_cols(dat_cis, common_cols)
gw_sub  <- if (nrow(gw_sub) > 0) ensure_cols(gw_sub, common_cols) else
             dat_cis[0]  # empty with same columns

dat_all <- rbind(gw_sub, dat_cis)
message("[DEBUG] dat_all: ", nrow(dat_all), " rows (", nrow(gw_sub), " GW + ",
        nrow(dat_cis), " cis)")

if (nrow(dat_all) == 0) {
  message("[warn] no valid MR results for pathway: ", pathway, " - saving placeholder")
  err_msgs <- if (nrow(dat_cis_err) > 0) {
    paste(dat_cis_err[, paste0(gwas_name, "/", gene_name, ": ", error)], collapse = "\n")
  } else { "No results found" }
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("No valid MR results\n(", pathway, ")\n\n", err_msgs),
             hjust = 0.5, vjust = 0.5, size = 4, colour = "grey40") +
    theme_void()
  ggsave(forest_plot, p_empty, width = 8, height = 4, dpi = 150, bg = "white")
  message("[DEBUG] placeholder saved, exiting cleanly")
  quit(save = "no", status = 0)
}

# OR cols
dat_all[, or   := exp(b)]
dat_all[, lo95 := exp(b - 1.96 * b_se)]
dat_all[, hi95 := exp(b + 1.96 * b_se)]
dat_all[, sig  := !is.na(p) & p < 0.05]

dat_all <- dat_all[outcome %in% outcome_order & method %in% method_order]
dat_all[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]
dat_all[, method_label  := factor(method_labels[method],   levels = method_labels)]

# ordering
dat_all[, exp_key := paste0(trimws(gwas_name), "__", trimws(gene_name))]
exp_order  <- unique(dat_all$exp_key)
n_exp      <- length(exp_order)
n_gw       <- sum(exp_order %in% paste0(unique(gw_sub$gwas_name), "__Genome-wide"))
message("[DEBUG] exposures (", n_exp, "): ", paste(exp_order, collapse = ", "))
message("[DEBUG] GW section: ", n_gw, " rows; cis section: ", n_exp - n_gw, " rows")

# plotting things
gap     <- 1.5
n_meth  <- length(method_order)
offsets <- setNames(
  seq(-(n_meth - 1) / 2, (n_meth - 1) / 2, by = 1) * 0.27,
  method_order
)

dat_all[, exp_idx := match(exp_key, exp_order)]
dat_all[, y_pos   := -(exp_idx * gap) + offsets[method]]

label_src <- dat_all[method == "mr_ivw" & outcome == "heart_failure",
                     .(gwas_name = first(gwas_name),
                       gene_name = first(gene_name),
                       n_snp     = first(n_snp),
                       fstat     = first(fstat)),
                     by = exp_key]
label_src_fb <- dat_all[method == "mr_ivw",
                        .(gwas_name = first(gwas_name),
                          gene_name = first(gene_name),
                          n_snp     = first(n_snp),
                          fstat     = first(fstat)),
                        by = exp_key]
label_src <- rbind(label_src, label_src_fb[!exp_key %in% label_src$exp_key])

label_src[, axis_label := paste0(
  fmt_name(gwas_name), " (", toupper(trimws(gene_name)), ")\n",
  "nSNP=", fifelse(is.na(n_snp), "?", as.character(n_snp)),
  "  F=",  fifelse(is.na(fstat), "?", as.character(round(fstat, 1)))
)]
label_src[, y_center := -(match(exp_key, exp_order) * gap)]

y_breaks <- label_src[match(exp_order, exp_key), y_center]
y_labels <- label_src[match(exp_order, exp_key), axis_label]

#  exposure band shading
shade_idx <- which(seq_along(exp_order) %% 2 == 0)
shade_dat <- if (length(shade_idx) > 0) {
  data.table(
    ymin = -(shade_idx * gap + gap / 2),
    ymax = -(shade_idx * gap - gap / 2)
  )
} else {
  data.table(ymin = numeric(0), ymax = numeric(0))
}

sep_y <- if (n_gw > 0 && n_gw < n_exp) -(n_gw * gap + gap / 2) else NULL

x_lo_cap  <- 0.05
x_hi_hard <- 3.0

caps <- dat_all[, {
  any_exceed <- any(!is.na(hi95) & hi95 > x_hi_hard)
  x_hi <- if (any_exceed) x_hi_hard else {
    max_val <- max(c(hi95, or), na.rm = TRUE)
    ceiling(max_val * 10) / 10 + 0.1
  }
  .(x_hi_cap = x_hi)
}, by = outcome]
message("[DEBUG] per-outcome x_hi_cap: ",
        paste(caps$outcome, round(caps$x_hi_cap, 2), sep = "=", collapse = ", "))

dat_all <- merge(dat_all, caps, by = "outcome")
dat_all[, lo95_clip := pmax(lo95, x_lo_cap)]
dat_all[, hi95_clip := pmin(hi95, x_hi_cap)]
dat_all[, trunc_lo  := !is.na(lo95) & lo95 < x_lo_cap]
dat_all[, trunc_hi  := !is.na(hi95) & hi95 > x_hi_cap]

# plot
plot_dat <- dat_all[!is.na(b) & !is.na(b_se)]
message("[DEBUG] plot_dat: ", nrow(plot_dat), " rows")
message("[DEBUG] truncated CIs: lo=", sum(plot_dat$trunc_lo, na.rm = TRUE),
        " hi=", sum(plot_dat$trunc_hi, na.rm = TRUE))
message("[DEBUG] building ggplot object...")

p <- ggplot(plot_dat,
            aes(x      = or,
                y      = y_pos,
                colour = method_label,
                shape  = method_label,
                size   = method_label)) +
  geom_rect(
    data        = shade_dat,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE,
    fill        = "grey96",
    colour      = NA
  ) +
  {if (!is.null(sep_y))
    geom_hline(yintercept = sep_y, colour = "grey60",
               linewidth = 0.6, linetype = "dashed")
  } +
  geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.45) +
  geom_errorbarh(
    aes(xmin = lo95_clip, xmax = hi95_clip),
    height    = 0.13,
    linewidth = 0.45,
    alpha     = 0.75
  ) +
  geom_segment(
    data        = plot_dat[trunc_hi == TRUE],
    aes(x       = x_hi_cap * 0.97,
        xend    = x_hi_cap,
        y       = y_pos,
        yend    = y_pos,
        colour  = method_label),
    arrow       = arrow(length = unit(0.12, "cm"), type = "closed"),
    linewidth   = 0.5,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data        = plot_dat[trunc_lo == TRUE],
    aes(x       = x_lo_cap * 1.5,
        xend    = x_lo_cap,
        y       = y_pos,
        yend    = y_pos,
        colour  = method_label),
    arrow       = arrow(length = unit(0.12, "cm"), type = "closed"),
    linewidth   = 0.5,
    inherit.aes = FALSE
  ) +
  geom_point(
    data  = plot_dat[or >= x_lo_cap & or <= x_hi_cap],
    alpha = 0.9
  ) +
  geom_point(
    data        = plot_dat[method == "mr_ivw" & sig == TRUE &
                           or >= x_lo_cap & or <= x_hi_cap],
    shape       = 21,
    size        = 4.5,
    stroke      = 0.7,
    colour      = "black",
    fill        = NA,
    inherit.aes = FALSE,
    aes(x = or, y = y_pos)
  ) +
  scale_colour_manual(name = "Method", values = method_colours) +
  scale_shape_manual( name = "Method", values = method_shapes)  +
  scale_size_manual(  name = "Method", values = method_sizes)   +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(add = gap * 0.6)
  ) +
  facet_wrap(~ outcome_label, nrow = 1, scales = "free_x") +
  labs(
    title = tools::toTitleCase(gsub("_", " ", pathway)),
    x     = "OR (95% CI)",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(colour = "grey70"),
    axis.text.y        = element_text(size = 8.5, lineheight = 0.85, hjust = 1),
    axis.ticks.y       = element_blank(),
    legend.position    = "bottom",
    legend.key.size    = unit(0.9, "lines"),
    legend.text        = element_text(size = 9),
    legend.box         = "horizontal",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin        = margin(8, 12, 8, 8)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3.5), nrow = 1),
    shape  = guide_legend(nrow = 1),
    size   = "none"
  )


fig_h <- max(4.0, n_exp * gap * 0.6 + 2.5)
fig_w <- 11.0
message("[DEBUG] saving to: ", forest_plot,
        " (", round(fig_w, 1), " x ", round(fig_h, 1), " in)")
ggsave(forest_plot, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
message("[DEBUG] done")
