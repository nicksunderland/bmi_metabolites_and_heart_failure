

# testing
if (FALSE) {
  model_dt_bbs_ms  <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "tmp_objects", "model_df_raw_bbs_metabolon.tsv")
  model_dt_bbs_nmr <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "tmp_objects", "model_df_raw_bbs_nightingale.tsv")
    # file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "tmp_objects", "model_df_raw_direct_metabolon.tsv")
    # file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "tmp_objects", "model_df_raw_direct_nightingale.tsv")
  metab_name_map <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "scripts", "gwas_metab_name_map.xlsx")
  out_dir <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "figures", "raw_metab_changes")
}


# requirements ----
library(data.table)
library(ggplot2)
library(readxl)
library(fs)
library(patchwork)


# read ----
bbs_ms  <- fread(model_dt_bbs_ms)
bbs_nmr <- fread(model_dt_bbs_nmr)
map     <- read_xlsx(metab_name_map, sheet = 1) |> as.data.table()


# clean ----
dat <- merge(bbs_ms [, .SD, .SDcols = c("study_id", "timepoint", "bmi", grep("^metab_", names(bbs_ms), value = TRUE))],
             bbs_nmr[, .SD, .SDcols = c("study_id", "timepoint", "bmi", grep("^metab_", names(bbs_nmr), value = TRUE))],
             by = c("study_id", "timepoint"), all = TRUE, suffixes = c("_ms", "_nmr"))
dat[, bmi := fcoalesce(bmi_ms, bmi_nmr)]
setcolorder(dat, c("study_id", "timepoint", "bmi"))
dat[, c("bmi_ms", "bmi_nmr") := NULL]
dat[, timepoint := factor(timepoint, levels = c("baseline", "end"))]




metabs <- grep("^metab_", names(dat), value = TRUE)


library(data.table)
library(ggplot2)
library(ggExtra)
library(patchwork)

for (i in seq_along(metabs)) {

  ## ---- metadata ----
  label <- map[id == metabs[i], label]
  pathway_group <- map[id == metabs[i], pathway_group]

  ## ---- long data ----
  metab_dat <- dat[, .(
    study_id,
    timepoint,
    bmi,
    metab = get(metabs[i])
  )]

  metab_dat[, timepoint := factor(timepoint, levels = c("baseline", "end"))]

  ## ---- compute deltas once ----
  delta_dat <- metab_dat[
    ,
    .(
      delta_metab = metab[timepoint == "end"] -
        metab[timepoint == "baseline"],
      delta_bmi   = bmi[timepoint == "end"] -
        bmi[timepoint == "baseline"]
    ),
    by = study_id
  ]

  ## attach direction back to long data
  metab_dat <- merge(
    metab_dat,
    delta_dat[, .(study_id, delta_metab, delta_bmi)],
    by = "study_id",
    all.x = TRUE
  )

  metab_dat[, direction := fifelse(delta_metab >= 0, "up", "down")]

  ## ===============================
  ## Plot 1: paired box + trajectories
  ## ===============================

  p1 <- ggplot(metab_dat, aes(x = timepoint, y = metab)) +

    # trajectories
    geom_line(
      aes(group = study_id, color = direction),
      alpha = 0.15,
      linewidth = 0.6
    ) +

    geom_point(
      aes(group = study_id),
      color = "grey40",
      alpha = 0.15,
      size = 1.5
    ) +

    stat_boxplot(
      geom = "errorbar",
      width = 0.25,
      linewidth = 0.6
    ) +

    # boxplots behind
    geom_boxplot(
      outliers = FALSE,
      fill = "white",
      alpha = 0.8,
      width = 0.5
    ) +



    scale_color_manual(
      values = c("up" = "darkgreen", "down" = "red3")
    ) +

    guides(color = "none") +

    labs(
      x = NULL,
      y = label
    ) +

    theme_classic()

  ## =================================
  ## Plot 2: Δmetabolite vs ΔBMI
  ## =================================

  p2 <- ggplot(
    delta_dat[!is.na(delta_metab) & abs(delta_bmi) > 0.5],
    aes(x = delta_bmi, y = delta_metab)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(color = "royalblue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
    labs(
      x = "\u0394 BMI",
      y = paste0("\u0394 ", label)
    ) +
    theme_classic()

  p2m <- ggMarginal(
    p2,
    type = "histogram",
    margins = "both",
    bins = 30,
    fill = "grey70",
    color = "grey30"
  )

  ## ===============================
  ## Combine
  ## ===============================

  comb <- p1 | p2m

  print(comb)

  ## ===============================
  ## Save
  ## ===============================
  fp_out <- file.path(
    out_dir,
    fs::path_sanitize(pathway_group),
    paste0(fs::path_sanitize(label), ".png")
  )

  dir.create(dirname(fp_out), showWarnings = FALSE, recursive = TRUE)

  ggsave(
    filename = fp_out,
    plot = comb,
    width = 9,
    height = 5,
    dpi = 300
  )
}





