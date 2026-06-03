#### DESCRIPTION #######################################
# Two-sample MR with unique instruments only.
# After harmonisation, removes any SNP that appears in
# more than one metabolite instrument (per snp_share_counts.tsv).
# Methods: Wald ratio and IVW only (no Steiger, no PRESSO).
########################################################

#### SET INPUT #########################################
instrument_file <- snakemake@input[["instrument_file"]]
outcome_file    <- snakemake@input[["outcome_file"]]
mapping_file    <- snakemake@input[["mapping_file"]]
snp_share_file  <- snakemake@input[["snp_share_file"]]
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
# output
results_file    <- snakemake@output[["results_file"]]
########################################################

if (FALSE) {
  # ── TESTING OVERRIDES ────────────────────────────────────────────────────────
  repo_dir        <- Sys.getenv("HF_METABOLITE_REPO2")
  instrument_file <- file.path(repo_dir, "output/tables/instruments/genome_wide/GCST90199629_instrument.tsv")
  outcome_file    <- Sys.getenv("HF_GWAS_PATH")
  mapping_file    <- file.path(repo_dir, "output/tables/global_variant_map.tsv.gz")
  snp_share_file  <- file.path(repo_dir, "output/tables/instruments/snp_share_counts.tsv")
  exp_id          <- "GCST90199629"
  exp_binary      <- FALSE
  exp_ncase       <- NULL
  exp_ncontrol    <- NULL
  exp_prevalence  <- NULL
  out_id          <- "heart_failure"
  out_build       <- "Hg19"
  out_binary      <- TRUE
  out_ncase       <- NULL
  out_ncontrol    <- NULL
  out_prevalence  <- NULL
  cores           <- 2L
  results_file    <- "/tmp/test_mr_unique_results.tsv"
  outcome_map <- list(
    chr  = list(name = "chr",  alias = "chromosome",            type = "character"),
    bp   = list(name = "bp",   alias = "base_pair_location",    type = "integer"),
    ea   = list(name = "ea",   alias = "effect_allele",         type = "character"),
    oa   = list(name = "oa",   alias = "other_allele",          type = "character"),
    beta = list(name = "beta", alias = "beta",                  type = "numeric"),
    se   = list(name = "se",   alias = "standard_error",        type = "numeric"),
    p    = list(name = "p",    alias = "p_value",               type = "numeric"),
    eaf  = list(name = "eaf",  alias = "effect_allele_frequency", type = "numeric")
  )
  # ─────────────────────────────────────────────────────────────────────────────
}

cat("=== INPUT SUMMARY ===\n")
cat("instrument_file:", instrument_file, "\n")
cat("outcome_file   :", outcome_file, "\n")
cat("mapping_file   :", mapping_file, "\n")
cat("snp_share_file :", snp_share_file, "\n")
cat("exp_id         :", exp_id, "\n")
cat("exp_binary     :", exp_binary, "\n")
cat("out_id         :", out_id, "\n")
cat("out_build      :", out_build, "\n")
cat("out_binary     :", out_binary, "\n")
cat("cores          :", cores, "\n")
cat("results_file   :", results_file, "\n")
cat("=====================\n")


library(data.table)
library(fst)
library(yaml)

setDTthreads(cores)

get_col <- function(map, key) map[[key]][["alias"]]

na_methods <- c("Wald ratio", "Inverse variance weighted")
fail_dt <- data.table(
  exp_id    = exp_id,
  out_id    = out_id,
  method    = na_methods,
  nsnp      = 0L,
  nsnp_removed_shared = NA_integer_,
  fstat     = NA_real_,
  b         = NA_real_,
  b_se      = NA_real_,
  p         = NA_real_
)

write_fail <- function(dt, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, file, sep = "\t")
}

# ── 1. READ INSTRUMENT ────────────────────────────────────────────────────────
cat("Reading instrument file:", basename(instrument_file), "\n")
exp <- fread(instrument_file)
cat(nrow(exp), "instruments\n")

# ── 2. READ SNP SHARE COUNTS ─────────────────────────────────────────────────
cat("Reading SNP share counts\n")
snp_share <- fread(snp_share_file)
shared_snps <- snp_share[n_metabolites > 1, rsid]
cat(length(shared_snps), "SNPs appear in >1 metabolite instrument\n")

# ── 3. READ VARIANT MAP ───────────────────────────────────────────────────────
cat("Reading variant map\n")
vmap <- fread(mapping_file, select = c("chr", "bp", "bp_b38", "ea", "oa", "rsid"))
vmap[, chr    := as.character(chr)]
vmap[, bp     := as.integer(bp)]
vmap[, bp_b38 := as.integer(bp_b38)]
vmap[, ea     := tolower(as.character(ea))]
vmap[, oa     := tolower(as.character(oa))]
vmap[, rsid   := as.character(rsid)]

# ── 4. READ OUTCOME GWAS ─────────────────────────────────────────────────────
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

old <- c(get_col(outcome_map, "chr"), get_col(outcome_map, "bp"),
         get_col(outcome_map, "ea"),  get_col(outcome_map, "oa"),
         get_col(outcome_map, "beta"),get_col(outcome_map, "se"),
         get_col(outcome_map, "p"),   get_col(outcome_map, "eaf"))
new <- c("chr", "bp", "ea", "oa", "beta", "se", "p", "eaf")
setnames(out, old, new)
if (!is.null(n_col) && n_col != "n") setnames(out, n_col, "n")

out[, chr  := as.character(chr)]
out[, bp   := as.integer(bp)]
out[, ea   := tolower(as.character(ea))]
out[, oa   := tolower(as.character(oa))]
out[, beta := as.numeric(beta)]
out[, se   := as.numeric(se)]
out[, p    := as.numeric(p)]
out[, eaf  := as.numeric(eaf)]

if (get_col(outcome_map, "p") == "neg_log_10_p_value") {
  cat("Converting neg_log10 p-value\n")
  out[, p := 10^(-p)]
}

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

n_before <- nrow(out)
out <- out[
  !is.na(chr) & !is.na(bp) & !is.na(ea) & !is.na(oa) &
    !is.na(beta) & !is.infinite(beta) & abs(beta) < 20 &
    !is.na(se)   & !is.infinite(se) & se > 0 &
    !is.na(p)    & !is.infinite(p) & p > 0 & p <= 1 &
    !is.na(eaf)  & eaf > 0 & eaf < 1
]
cat("Filtered", n_before - nrow(out), "rows,", nrow(out), "remaining in outcome\n")

# ── 5. ANNOTATE OUTCOME WITH b37 RSID ────────────────────────────────────────
cat("Annotating outcome with b37 rsid\n")
if (out_build == "Hg38") {
  out[, bp_b38 := bp]
  out[vmap, c("rsid", "bp") := list(i.rsid, i.bp), on = c("chr", "bp_b38", "ea", "oa")]
  out[, bp_b38 := NULL]
} else {
  out[vmap, rsid := i.rsid, on = c("chr", "bp", "ea", "oa")]
}
out[is.na(rsid) | !grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]
cat("Outcome RSIDs annotated:", out[grepl("^rs[0-9]+$", rsid), .N], "/", nrow(out), "\n")

# ── 6. SUBSET OUTCOME TO INSTRUMENT SNPS ─────────────────────────────────────
cat("Subsetting outcome to instrument SNPs\n")
out <- out[rsid %in% exp$rsid]
cat(nrow(out), "outcome SNPs matching instruments\n")

if (nrow(out) == 0) {
  cat("WARNING: No instrument SNPs found in outcome GWAS for", exp_id, "->", out_id, "\n")
  write_fail(fail_dt, results_file)
  quit(save = "no", status = 0)
}

# ── 7. FORMAT AND HARMONISE ───────────────────────────────────────────────────
cat("Formatting and harmonising via TwoSampleMR\n")
library(TwoSampleMR)

exposure_dat <- TwoSampleMR::format_data(
  dat               = as.data.frame(exp),
  type              = "exposure",
  snp_col           = "rsid",
  beta_col          = "beta",
  se_col            = "se",
  pval_col          = "p",
  eaf_col           = "eaf",
  effect_allele_col = "ea",
  other_allele_col  = "oa",
  chr_col           = "chr",
  pos_col           = "bp",
  samplesize_col    = "n"
)
exposure_dat$id.exposure <- exp_id

outcome_dat <- TwoSampleMR::format_data(
  dat               = as.data.frame(out),
  type              = "outcome",
  snp_col           = "rsid",
  beta_col          = "beta",
  se_col            = "se",
  pval_col          = "p",
  eaf_col           = "eaf",
  effect_allele_col = "ea",
  other_allele_col  = "oa",
  chr_col           = "chr",
  pos_col           = "bp",
  samplesize_col    = "n"
)
outcome_dat$id.outcome <- out_id

har_dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2)
har_dat <- har_dat[har_dat$mr_keep, ]
cat(nrow(har_dat), "SNPs retained after harmonisation\n")

if (nrow(har_dat) == 0) {
  cat("WARNING: No SNPs retained after harmonisation for", exp_id, "->", out_id, "\n")
  write_fail(fail_dt, results_file)
  quit(save = "no", status = 0)
}

# ── 8. REMOVE SHARED SNPs ────────────────────────────────────────────────────
n_before_filter <- nrow(har_dat)
har_dat <- har_dat[!har_dat$SNP %in% shared_snps, ]
nsnp_removed <- n_before_filter - nrow(har_dat)
cat(nsnp_removed, "shared SNPs removed,", nrow(har_dat), "unique SNPs remaining\n")

if (nrow(har_dat) == 0) {
  cat("WARNING: No unique SNPs remaining after shared SNP removal for", exp_id, "->", out_id, "\n")
  fail_dt[, nsnp_removed_shared := nsnp_removed]
  write_fail(fail_dt, results_file)
  quit(save = "no", status = 0)
}

fstat <- sum((har_dat$beta.exposure / har_dat$se.exposure)^2) / nrow(har_dat)
cat("Mean F-statistic:", round(fstat, 2), "\n")

# ── 9. RUN MR (Wald ratio / IVW only) ────────────────────────────────────────
cat("Running MR methods\n")
res_mr <- TwoSampleMR::mr(har_dat, method_list = c("mr_wald_ratio", "mr_ivw"))

res_main <- as.data.table(res_mr)[, .(method, nsnp, b, se, pval)]
setnames(res_main, c("se", "pval"), c("b_se", "p"))

res_main[, exp_id             := exp_id]
res_main[, out_id             := out_id]
res_main[, fstat              := fstat]
res_main[, nsnp_removed_shared := nsnp_removed]
setcolorder(res_main, c("exp_id", "out_id", "method", "nsnp", "nsnp_removed_shared", "fstat"))

cat("Writing:", results_file, "\n")
dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
fwrite(res_main, results_file, sep = "\t")
