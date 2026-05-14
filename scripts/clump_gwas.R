#### DESCRIPTION #######################################
# This script clumps a GWAS to create a genome-wide
# instrument, a set of independent variants for a given trait.

#### SET INPUT #########################################
gwas_file    <- snakemake@input[["gwas_file"]]
mapping_file <- snakemake@input[["mapping_file"]]
map          <- snakemake@params[["map"]]
build        <- snakemake@params[["build"]]
clump_r2     <- snakemake@params[["r2"]]
clump_p1     <- snakemake@params[["p1"]]
clump_kb     <- snakemake@params[["kb"]]
pfile        <- snakemake@params[["pfile"]]
trait        <- snakemake@params[["trait"]]
id           <- snakemake@params[["id"]]
cores        <- snakemake@resources[["cpus_per_task"]]
clump_file   <- snakemake@output[["clump_file"]]
########################################################

cat("=== INPUT SUMMARY ===\n")
cat("gwas_file   :", gwas_file, "\n")
cat("mapping_file:", mapping_file, "\n")
cat("build       :", build, "\n")
cat("trait       :", trait, "\n")
cat("id          :", id, "\n")
cat("pfile       :", pfile, "\n")
cat("clump_r2    :", clump_r2, "\n")
cat("clump_p1    :", clump_p1, "\n")
cat("clump_kb    :", clump_kb, "\n")
cat("=====================\n")

library(data.table)
library(fst)
library(yaml)

# set threads for data.table
setDTthreads(cores)
cat("Threads:", getDTthreads(), "\n")

# helper to get alias from map
get_col <- function(map, key) map[[key]][["alias"]]

# read GWAS file
cat("Reading GWAS file:", basename(gwas_file), "\n")
cols <- c(
  get_col(map, "chr"),
  get_col(map, "bp"),
  get_col(map, "ea"),
  get_col(map, "oa"),
  get_col(map, "beta"),
  get_col(map, "se"),
  get_col(map, "p"),
  get_col(map, "eaf")
)
# add n if in map
n_col <- if (!is.null(map[["n"]])) get_col(map, "n") else NULL
if (!is.null(n_col)) cols <- c(cols, n_col)

if (grepl("\\.fst$", gwas_file)) {
  gwas <- fst::read_fst(gwas_file, columns = cols, as.data.table = TRUE)
} else {
  gwas <- data.table::fread(gwas_file, select = cols)
}

# rename to standard names
std_names <- c("chr", "bp", "ea", "oa", "beta", "se", "p", "eaf")
old_names <- c(
  get_col(map, "chr"), get_col(map, "bp"),
  get_col(map, "ea"),  get_col(map, "oa"),
  get_col(map, "beta"),get_col(map, "se"),
  get_col(map, "p"),   get_col(map, "eaf")
)
setnames(gwas, old_names, std_names)
if (!is.null(n_col) && n_col != "n") setnames(gwas, n_col, "n")

# convert neg_log10 p-value if needed
if (get_col(map, "p") == "neg_log_10_p_value") {
  cat("Converting neg_log10 p-value to p-value\n")
  gwas[, p := 10^(-p)]
}

# ensure types
gwas[, chr  := as.character(chr)]
gwas[, bp   := as.integer(bp)]
gwas[, ea   := toupper(as.character(ea))]
gwas[, oa   := toupper(as.character(oa))]
gwas[, beta := as.numeric(beta)]
gwas[, se   := as.numeric(se)]
gwas[, p    := as.numeric(p)]
gwas[, eaf  := as.numeric(eaf)]

# get n from yaml if not in gwas
if (!"n" %in% names(gwas)) {
  cat("N not in map, trying yaml metadata\n")
  yaml_file <- sub("\\.fst$|\\.tsv\\.gz$|\\.tsv$", ".tsv-meta.yaml", gwas_file)
  if (file.exists(yaml_file)) {
    meta    <- yaml::read_yaml(yaml_file)
    n_total <- tryCatch({
      sizes <- sapply(meta[["samples"]], function(s) s[["sample_size"]])
      sizes <- sizes[!sapply(sizes, is.null)]
      if (length(sizes) == 0) NA_integer_ else as.integer(sum(unlist(sizes)))
    }, error = function(e) {
      cat("Could not parse sample_size from yaml:", conditionMessage(e), "\n")
      NA_integer_
    })
    cat("N from yaml:", n_total, "\n")
    gwas[, n := n_total]
  } else {
    cat("No yaml found, N will be NA\n")
    gwas[, n := NA_integer_]
  }
}

# basic filters
n_before <- nrow(gwas)
gwas <- gwas[
  !is.na(chr) & !is.na(bp) & !is.na(ea) & !is.na(oa) &
  !is.na(beta) & !is.infinite(beta) & abs(beta) < 20 &
  !is.na(se) & !is.infinite(se) &
  !is.na(p) & !is.infinite(p) &
  !is.na(eaf) & eaf > 0 & eaf < 1
]
cat("Filtered", n_before - nrow(gwas), "rows,", nrow(gwas), "remaining\n")

# join variant map - always use b37 positions and rsid for clumping
cat("Reading variant map:", mapping_file, "\n")
vmap <- fread(mapping_file, select = c("rsid", "chr", "bp", "bp_b38", "ea", "oa"))

# join on native build coordinates to get b37 position and rsid
if (build == "Hg38") {
  gwas[, bp_b38 := bp]
  gwas[vmap, c("rsid", "bp") := list(i.rsid, i.bp), on = c("chr", "bp_b38", "ea", "oa")]
  gwas[, bp_b38 := NULL]
} else {
  gwas[vmap, rsid := i.rsid, on = c("chr", "bp", "ea", "oa")]
}

# fallback rsid
gwas[is.na(rsid) | !grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]
cat("RSIDs annotated:", gwas[grepl("^rs[0-9]+$", rsid), .N], "/", nrow(gwas), "\n")

# clump
cat("Clumping\n")
clumps <- genepi.utils::clump(gwas      = gwas,
                              p1        = clump_p1,
                              r2        = clump_r2,
                              kb        = clump_kb,
                              plink2    = NULL,
                              plink_ref = pfile)[index == TRUE, ]

# clean up and write
clumps[, c("index", "clump", "trait", "id") := list(NULL, NULL, trait, id)]
cat("Writing:", clump_file, "\n")
data.table::fwrite(clumps, clump_file, sep = "\t")