#### DESCRIPTION #######################################
# Generalised two-sample MR script
# Exposure: pre-clumped instrument file (b37)
# Outcome:  any GWAS file with column map
# Build:    always harmonised to b37 via variant map

#### SET INPUT #########################################
instrument_file <- snakemake@input[["instrument_file"]]
outcome_file    <- snakemake@input[["outcome_file"]]
mapping_file    <- snakemake@input[["mapping_file"]]
# exposure params
exp_id          <- snakemake@params[["exp_id"]]
exp_binary      <- snakemake@params[["exp_binary"]]
exp_ncase       <- snakemake@params[["exp_ncase"]]
exp_ncontrol    <- snakemake@params[["exp_ncontrol"]]
exp_prevalence  <- snakemake@params[["exp_prevalence"]]
# outcome params
out_id          <- snakemake@params[["out_id"]]
outcome_map     <- snakemake@params[["outcome_map"]]
out_build       <- snakemake@params[["out_build"]]
out_binary      <- snakemake@params[["out_binary"]]
out_ncase       <- snakemake@params[["out_ncase"]]
out_ncontrol    <- snakemake@params[["out_ncontrol"]]
out_prevalence  <- snakemake@params[["out_prevalence"]]
# other params
cores           <- snakemake@resources[["cpus_per_task"]]
steiger_filter  <- snakemake@params[["steiger_filter"]]
# output
results_file    <- snakemake@output[["results_file"]]
########################################################

if (FALSE) {
# ── TESTING OVERRIDES (comment out when running via snakemake) ───────────────
instrument_file <- "/Users/xx20081/git/bmi_metabolites_and_heart_failure/output/tables/instruments/genome_wide/bmi_instrument.tsv"
outcome_file    <- "/Users/xx20081/Downloads/GCST90200363_buildGRCh38.tsv.gz"
mapping_file    <- "/Users/xx20081/git/bmi_metabolites_and_heart_failure/output/tables/global_variant_map.tsv.gz"
exp_id          <- "bmi"
exp_binary      <- FALSE
exp_ncase       <- NULL
exp_ncontrol    <- NULL
exp_prevalence  <- NULL
out_id          <- "GCST90200363"
out_build       <- "Hg38"
out_binary      <- FALSE
out_ncase       <- NULL
out_ncontrol    <- NULL
out_prevalence  <- NULL
steiger_filter  <- TRUE
cores           <- 4L
results_file    <- "/tmp/test_mr_results.tsv"
outcome_map <- list(
  chr  = list(name = "chr",  alias = "chromosome",  type = "character"),
  bp   = list(name = "bp",   alias = "base_pair_location",  type = "integer"),
  ea   = list(name = "ea",   alias = "effect_allele",   type = "character"),
  oa   = list(name = "oa",   alias = "other_allele",   type = "character"),
  beta = list(name = "beta", alias = "beta", type = "numeric"),
  se   = list(name = "se",   alias = "standard_error",   type = "numeric"),
  p    = list(name = "p",    alias = "p_value",    type = "numeric"),
  eaf  = list(name = "eaf",  alias = "effect_allele_frequency",  type = "numeric")
)
# ─────────────────────────────────────────────────────────────────────────────
}

cat("=== INPUT SUMMARY ===\n")
cat("instrument_file:", instrument_file, "\n")
cat("outcome_file   :", outcome_file, "\n")
cat("mapping_file   :", mapping_file, "\n")
cat("exp_id         :", exp_id, "\n")
cat("exp_binary     :", exp_binary, "\n")
cat("exp_ncase      :", exp_ncase, "\n")
cat("exp_ncontrol   :", exp_ncontrol, "\n")
cat("exp_prevalence :", exp_prevalence, "\n")
cat("out_id         :", out_id, "\n")
cat("out_build      :", out_build, "\n")
cat("out_binary     :", out_binary, "\n")
cat("out_ncase      :", out_ncase, "\n")
cat("out_ncontrol   :", out_ncontrol, "\n")
cat("out_prevalence :", out_prevalence, "\n")
cat("steiger_filter :", steiger_filter, "\n")
cat("cores          :", cores, "\n")
cat("results_file   :", results_file, "\n")
cat("=====================\n")


library(data.table)
library(fst)
library(yaml)

setDTthreads(cores)
cat("Threads:", getDTthreads(), "\n")

get_col <- function(map, key) map[[key]][["alias"]]

# ── Canonical NA result — updated in-place before each early exit ─────────────
na_methods <- c(
  "Inverse variance weighted", "MR Egger", "Weighted median",
  "Weighted mode", "MR-PRESSO-raw", "MR-PRESSO-corrected"
)
fail_dt <- data.table(
  exp_id                = exp_id,
  out_id                = out_id,
  method                = na_methods,
  nsnp                  = 0L,
  fstat                 = NA_real_,
  b                     = NA_real_,
  b_se                  = NA_real_,
  p                     = NA_real_,
  Q                     = NA_real_,
  Q_df                  = NA_integer_,
  Q_pval                = NA_real_,
  egger_intercept       = NA_real_,
  egger_intercept_se    = NA_real_,
  egger_intercept_p     = NA_real_,
  presso_global_p       = NA_character_,
  steiger_filtering     = steiger_filter,
  nsnp_steiger_filtered = NA_integer_,
  snps_steiger_filtered = NA_character_,
  exp_snp_r2            = NA_real_,
  out_snp_r2            = NA_real_,
  steiger_correct_dir   = NA,
  steiger_p             = NA_real_
)

# ── 1. READ INSTRUMENT (already b37, clumped) ────────────────────────────────
cat("Reading instrument file:", basename(instrument_file), "\n")
exp <- fread(instrument_file)
cat(nrow(exp), "instruments\n")

# ── 2. READ VARIANT MAP (for outcome rsid/b37 annotation) ────────────────────
cat("Reading variant map\n")
if (out_build == "Hg19") {
  vmap <- fread(mapping_file, select = c("chr", "bp", "bp_b38", "ea", "oa", "rsid"))
} else {
  vmap <- fread(mapping_file, select = c("chr", "bp", "bp_b38", "ea", "oa", "rsid"))
}
# type coercion
vmap[, chr       := as.character(chr)]
vmap[, bp        := as.integer(bp)]
vmap[, bp_b38    := as.integer(bp_b38)]
vmap[, ea        := tolower(as.character(ea))]
vmap[, oa        := tolower(as.character(oa))]
vmap[, rsid      := as.character(rsid)]

# ── 3. READ OUTCOME GWAS ─────────────────────────────────────────────────────
cat("Reading outcome file:", basename(outcome_file), "\n")
out_cols <- c(
  get_col(outcome_map, "chr"),
  get_col(outcome_map, "bp"),
  get_col(outcome_map, "ea"),
  get_col(outcome_map, "oa"),
  get_col(outcome_map, "beta"),
  get_col(outcome_map, "se"),
  get_col(outcome_map, "p"),
  get_col(outcome_map, "eaf")
)
n_col <- if (!is.null(outcome_map[["n"]])) get_col(outcome_map, "n") else NULL
if (!is.null(n_col)) out_cols <- c(out_cols, n_col)

if (grepl("\\.fst$", outcome_file)) {
  out <- read_fst(outcome_file, columns = out_cols, as.data.table = TRUE)
} else {
  out <- fread(outcome_file, select = out_cols)
}

# rename to standard
old <- c(get_col(outcome_map, "chr"), get_col(outcome_map, "bp"),
         get_col(outcome_map, "ea"),  get_col(outcome_map, "oa"),
         get_col(outcome_map, "beta"),get_col(outcome_map, "se"),
         get_col(outcome_map, "p"),   get_col(outcome_map, "eaf"))
new <- c("chr", "bp", "ea", "oa", "beta", "se", "p", "eaf")
setnames(out, old, new)
if (!is.null(n_col) && n_col != "n") setnames(out, n_col, "n")

# type coercion
out[, chr  := as.character(chr)]
out[, bp   := as.integer(bp)]
out[, ea   := tolower(as.character(ea))]
out[, oa   := tolower(as.character(oa))]
out[, beta := as.numeric(beta)]
out[, se   := as.numeric(se)]
out[, p    := as.numeric(p)]
out[, eaf  := as.numeric(eaf)]

# convert neg log10 p if needed
if (get_col(outcome_map, "p") == "neg_log_10_p_value") {
  cat("Converting neg_log10 p-value\n")
  out[, p := 10^(-p)]
}

# get n from yaml if missing
if (!"n" %in% names(out)) {
  yaml_file <- sub("\\.fst$|\\.tsv\\.gz$|\\.tsv$", ".tsv-meta.yaml", outcome_file)
  if (file.exists(yaml_file)) {
    meta    <- yaml::read_yaml(yaml_file)
    n_total <- tryCatch({
      sizes <- sapply(meta[["samples"]], function(s) s[["sample_size"]])
      sizes <- sizes[!sapply(sizes, is.null)]
      if (length(sizes) == 0) NA_integer_ else as.integer(sum(unlist(sizes)))
    }, error = function(e) NA_integer_)
    cat("N from yaml:", n_total, "\n")
    out[, n := n_total]
  } else {
    cat("No yaml found, N will be NA\n")
    out[, n := NA_integer_]
  }
}

# basic filters
n_before <- nrow(out)
out <- out[
  !is.na(chr) & !is.na(bp) & !is.na(ea) & !is.na(oa) &
  !is.na(beta) & !is.infinite(beta) & abs(beta) < 20 &
  !is.na(se)   & !is.infinite(se) & se > 0 &
  !is.na(p)    & !is.infinite(p) & p > 0 & p <= 1 &
  !is.na(eaf)  & eaf > 0 & eaf < 1
]
cat("Filtered", n_before - nrow(out), "rows,", nrow(out), "remaining in outcome\n")

# ── 4. ANNOTATE OUTCOME WITH b37 RSID AND POSITION ───────────────────────────
cat("Annotating outcome with b37 rsid and position\n")
if (out_build == "Hg38") {
  out[, bp_b38 := bp]
  out[vmap, c("rsid", "bp") := list(i.rsid, i.bp), on = c("chr", "bp_b38", "ea", "oa")]
  out[, bp_b38 := NULL]
} else {
  out[vmap, rsid := i.rsid, on = c("chr", "bp", "ea", "oa")]
}

# fallback rsid
out[is.na(rsid) | !grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]

cat("Outcome RSIDs annotated:", out[grepl("^rs[0-9]+$", rsid), .N], "/", nrow(out), "\n")


# ── 5. SUBSET OUTCOME TO INSTRUMENT SNPS ─────────────────────────────────────
cat("Subsetting outcome to instrument SNPs\n")
out <- out[rsid %in% exp$rsid]
cat(nrow(out), "outcome SNPs matching instruments\n")

if (nrow(out) == 0) {
  cat("WARNING: No instrument SNPs found in outcome GWAS for", exp_id, "->", out_id, "\n")
  fwrite(fail_dt, results_file, sep = "\t")
  quit(save = "no", status = 0)
}

# ── 6. FORMAT AND HARMONISE VIA TWOSAMPLEMR ──────────────────────────────────
cat("Formatting and harmonising via TwoSampleMR\n")
library(TwoSampleMR)

exposure_dat <- TwoSampleMR::format_data(
  dat                = as.data.frame(exp),
  type               = "exposure",
  snp_col            = "rsid",
  beta_col           = "beta",
  se_col             = "se",
  pval_col           = "p",
  eaf_col            = "eaf",
  effect_allele_col  = "ea",
  other_allele_col   = "oa",
  chr_col            = "chr",
  pos_col            = "bp",
  samplesize_col     = "n",
)
exposure_dat$id.exposure <- exp_id

outcome_dat <- TwoSampleMR::format_data(
  dat                = as.data.frame(out),
  type               = "outcome",
  snp_col            = "rsid",
  beta_col           = "beta",
  se_col             = "se",
  pval_col           = "p",
  eaf_col            = "eaf",
  effect_allele_col  = "ea",
  other_allele_col   = "oa",
  chr_col            = "chr",
  pos_col            = "bp",
  samplesize_col     = "n"
)
outcome_dat$id.outcome <- out_id

# harmonise - action=2 handles palindromic SNPs by EAF
har_dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2)
har_dat <- har_dat[har_dat$mr_keep, ]
cat(nrow(har_dat), "SNPs retained after harmonisation\n")

if (nrow(har_dat) == 0) {
  cat("WARNING: No SNPs retained after harmonisation for", exp_id, "->", out_id, "\n")
  fwrite(fail_dt, results_file, sep = "\t")
  quit(save = "no", status = 0)
}


# ── 7. STEIGER r², FILTERING, AND MR METHODS ─────────────────────────────────
cat("Calculating r² for Steiger test\n")

if (isTRUE(exp_binary)) {
  har_dat$units.exposure      <- "log odds"
  har_dat$ncase.exposure      <- exp_ncase
  har_dat$ncontrol.exposure   <- exp_ncontrol
  har_dat$prevalence.exposure <- exp_prevalence
  har_dat$r.exposure          <- TwoSampleMR::get_r_from_lor(
    lor        = har_dat$beta.exposure,
    af         = har_dat$eaf.exposure,
    ncase      = exp_ncase,
    ncontrol   = exp_ncontrol,
    prevalence = exp_prevalence
  )
} else {
  har_dat$units.exposure <- "SD"
  har_dat$r.exposure     <- TwoSampleMR::get_r_from_bsen(
    b  = har_dat$beta.exposure,
    se = har_dat$se.exposure,
    n  = har_dat$samplesize.exposure
  )
}

if (isTRUE(out_binary)) {
  har_dat$units.outcome      <- "log odds"
  har_dat$ncase.outcome      <- out_ncase
  har_dat$ncontrol.outcome   <- out_ncontrol
  har_dat$prevalence.outcome <- out_prevalence
  har_dat$r.outcome          <- TwoSampleMR::get_r_from_lor(
    lor        = har_dat$beta.outcome,
    af         = har_dat$eaf.outcome,
    ncase      = out_ncase,
    ncontrol   = out_ncontrol,
    prevalence = out_prevalence
  )
} else {
  har_dat$units.outcome <- "SD"
  har_dat$r.outcome     <- TwoSampleMR::get_r_from_bsen(
    b  = har_dat$beta.outcome,
    se = har_dat$se.outcome,
    n  = har_dat$samplesize.outcome
  )
}

har_dat <- TwoSampleMR::steiger_filtering(har_dat)

nsnp_steiger_filtered <- 0L
snps_steiger_filtered <- NA_character_
if (isTRUE(steiger_filter)) {
  steiger_fail_idxs <- which(har_dat$steiger_dir == FALSE & har_dat$steiger_pval < 0.05)
  nsnp_steiger_filtered <- length(steiger_fail_idxs)
  if (nsnp_steiger_filtered > 0) {
    snps_steiger_filtered <- paste(
      paste0(har_dat$SNP[steiger_fail_idxs], "_",
             har_dat$effect_allele.exposure[steiger_fail_idxs], "_",
             har_dat$other_allele.exposure[steiger_fail_idxs]),
      collapse = ";"
    )
    har_dat <- har_dat[-steiger_fail_idxs, ]
    cat(nsnp_steiger_filtered, "SNPs removed by Steiger filter,", nrow(har_dat), "remaining\n")
  }
}
fail_dt[, nsnp_steiger_filtered := nsnp_steiger_filtered]
fail_dt[, snps_steiger_filtered := snps_steiger_filtered]

if (nrow(har_dat) == 0) {
  cat("WARNING: No SNPs retained after Steiger filtering for", exp_id, "->", out_id, "\n")
  fwrite(fail_dt, results_file, sep = "\t")
  quit(save = "no", status = 0)
}

fstat <- sum((har_dat$beta.exposure / har_dat$se.exposure)^2) / nrow(har_dat)
cat("Mean F-statistic:", round(fstat, 2), "\n")

# ── 8. RUN MR METHODS ────────────────────────────────────────────────────────
cat("Running MR methods\n")

res_mr <- TwoSampleMR::mr(har_dat, method_list = c(
  "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"
))

res_het <- tryCatch({
  TwoSampleMR::mr_heterogeneity(har_dat, method_list = c("mr_ivw", "mr_egger_regression"))
}, error = function(e) {
  cat("Heterogeneity analysis failed:", conditionMessage(e), "\n"); NULL
})

res_pleio <- tryCatch({
  TwoSampleMR::mr_pleiotropy_test(har_dat)
}, error = function(e) {
  cat("Pleiotropy test failed:", conditionMessage(e), "\n"); NULL
})

res_presso <- tryCatch({
  TwoSampleMR::run_mr_presso(har_dat, NbDistribution = 1000, SignifThreshold = 0.05)[[1]]
}, error = function(e) {
  cat("MR-PRESSO failed:", conditionMessage(e), "\n"); NULL
})

# Overall Steiger directionality — run on post-filter instruments
res_steiger <- tryCatch({
  TwoSampleMR::directionality_test(har_dat)
}, error = function(e) {
  cat("Steiger directionality failed:", conditionMessage(e), "\n"); NULL
})


# ── 9. COMBINE AND WRITE ─────────────────────────────────────────────────────
cat("Combining results\n")

# main MR results
res_main <- as.data.table(res_mr)[, .(method, nsnp, b, se, pval)]
setnames(res_main, c("se", "pval"), c("b_se", "p"))

# heterogeneity - wide join onto main
if (!is.null(res_het) && is.data.frame(res_het) && nrow(res_het) > 0 && "method" %in% names(res_het)) {
  res_het_dt <- as.data.table(res_het)[, .(method, Q, Q_df, Q_pval)]
  res_main   <- merge(res_main, res_het_dt, by = "method", all.x = TRUE)
} else {
  cat("Heterogeneity results unavailable; setting Q columns to NA\n")
  res_main[, c("Q", "Q_df", "Q_pval") := list(NA_real_, NA_integer_, NA_real_)]
}

# egger intercept - broadcast to all rows
if (!is.null(res_pleio) && is.data.frame(res_pleio) && nrow(res_pleio) > 0) {
  res_main[, egger_intercept    := res_pleio$egger_intercept]
  res_main[, egger_intercept_se := res_pleio$se]
  res_main[, egger_intercept_p  := res_pleio$pval]
} else {
  res_main[, c("egger_intercept", "egger_intercept_se", "egger_intercept_p") := list(NA_real_, NA_real_, NA_real_)]
}


if (!is.null(res_presso)) {
  presso_dt <- tryCatch({

    raw           <- res_presso$`Main MR results`
    global        <- res_presso$`MR-PRESSO results`$`Global Test`
    outlier_obj   <- res_presso$`MR-PRESSO results`$`Distortion Test`
    n_outliers    <- if (!is.null(outlier_obj)) length(outlier_obj$`Outliers Indices`) else 0
    prop_outliers <- n_outliers / nrow(har_dat)
    distortion_p  <- if (!is.null(outlier_obj))  as.character(outlier_obj$Pvalue) else NA_character_

    r1 <- data.table(
      method = "MR-PRESSO-raw",
      nsnp   = sum(har_dat$mr_keep, na.rm = T),
      b      = as.numeric(raw[1, "Causal Estimate"]),
      b_se   = as.numeric(raw[1, "Sd"]),
      p      = as.numeric(raw[1, "P-value"]),
      presso_global_p = as.character(global$Pvalue),
      n_outliers      = n_outliers,
      prop_outliers   = prop_outliers,
      distortion_p    = NA_character_
    )
    r2 <- data.table(
      method = "MR-PRESSO-corrected",
      nsnp   = sum(har_dat$mr_keep, na.rm = T) - n_outliers,
      b      = as.numeric(raw[2, "Causal Estimate"]),
      b_se   = as.numeric(raw[2, "Sd"]),
      p      = as.numeric(raw[2, "P-value"]),
      presso_global_p = NA_character_,
      n_outliers      = NA_integer_,
      prop_outliers   = NA_real_,
      distortion_p    = distortion_p
    )
    rbind(r1, r2)

  }, error = function(e) {
    cat("PRESSO extraction failed:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(presso_dt)) {
    res_main <- rbindlist(list(res_main, presso_dt), fill = TRUE)
  }
}


# steiger
res_main[, steiger_filtering     := steiger_filter]
res_main[, nsnp_steiger_filtered := nsnp_steiger_filtered]
res_main[, snps_steiger_filtered := snps_steiger_filtered]
if (!is.null(res_steiger)) {
  res_main[, exp_snp_r2          := res_steiger$snp_r2.exposure]
  res_main[, out_snp_r2          := res_steiger$snp_r2.outcome]
  res_main[, steiger_correct_dir := res_steiger$correct_causal_direction]
  res_main[, steiger_p           := res_steiger$steiger_pval]
} else {
  res_main[, c("exp_snp_r2", "out_snp_r2", "steiger_correct_dir", "steiger_p") :=
             list(NA_real_, NA_real_, NA, NA_real_)]
}

# add ids and instrument strength
res_main[, exp_id := exp_id]
res_main[, out_id := out_id]
res_main[, fstat  := fstat]
setcolorder(res_main, c("exp_id", "out_id", "method", "nsnp", "fstat"))

cat("Writing:", results_file, "\n")
fwrite(res_main, results_file, sep = "\t")
