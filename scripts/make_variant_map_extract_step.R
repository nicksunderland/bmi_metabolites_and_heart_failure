#### DESCRIPTION #######################################
# Extracts novel variants (chr/bp/ea/oa) from a GWAS file
# by antijoining against the appropriate build base file

#### SET INPUT #########################################
input_gwas_fp    <- snakemake@input[["input_gwas_fp"]]
input_gwas_build <- snakemake@params[["input_gwas_build"]]
input_gwas_map   <- snakemake@params[["input_gwas_map"]]
base_hg19_fp     <- snakemake@params[["base_hg19_fp"]]
base_hg38_fp     <- snakemake@params[["base_hg38_fp"]]
base_hg19_map    <- snakemake@params[["base_hg19_map"]]
base_hg38_map    <- snakemake@params[["base_hg38_map"]]
cores            <- as.integer(snakemake@resources[["cpus_per_task"]])
out              <- snakemake@output[["out"]]
########################################################

# cat("=== INPUT SUMMARY ===\n")
# cat("input_gwas_fp   :", input_gwas_fp, "\n")
# cat("input_gwas_build:", input_gwas_build, "\n")
# cat("base_hg19_fp    :", base_hg19_fp, "\n")
# cat("base_hg38_fp    :", base_hg38_fp, "\n")
# cat("out             :", out, "\n")
# cat("input_gwas_map  :\n"); str(input_gwas_map)
# cat("base_hg19_map   :\n"); str(base_hg19_map)
# cat("base_hg38_map   :\n"); str(base_hg38_map)
# cat("=====================\n")


library(data.table)
library(fst)

setDTthreads(cores)
cat("Threads:", getDTthreads(), "\n")

# select correct base for this build
if (input_gwas_build == "Hg19") {
  base_fp  <- base_hg19_fp
  base_map <- base_hg19_map
} else {
  base_fp  <- base_hg38_fp
  base_map <- base_hg38_map
}

# standardised column names using map
get_col <- function(map, key) map[[key]][["alias"]]

# read input gwas - select and rename to standard names
cat("Reading GWAS file:", input_gwas_fp, "\n")
cols <- c(
  get_col(input_gwas_map, "chr"),
  get_col(input_gwas_map, "bp"),
  get_col(input_gwas_map, "ea"),
  get_col(input_gwas_map, "oa")
)

read_gwas <- function(fp, columns) {
  dt <- if (grepl("\\.fst$", fp)) {
    read_fst(fp, columns = columns, as.data.table = TRUE)
  } else {
    fread(fp, select = columns)
  }
  # complete cases only
  n_before <- nrow(dt)
  dt <- dt[complete.cases(dt)]
  n_dropped <- n_before - nrow(dt)
  if (n_dropped > 0) cat("Dropped", n_dropped, "rows with missing values from", fp, "\n")
  dt
}

gwas <- read_gwas(input_gwas_fp, cols)
setnames(gwas, cols, c("chr", "bp", "ea", "oa"))
gwas[, chr := as.character(chr)]
gwas[, bp  := as.integer(bp)]
gwas[, ea  := as.character(ea)]
gwas[, oa  := as.character(oa)]
gwas <- unique(gwas)


# antijoin against base if this is not the base file itself
if (input_gwas_fp != base_fp) {
  cat("Antijoining against base:", base_fp, "\n")

  base_cols <- c(
    get_col(base_map, "chr"),
    get_col(base_map, "bp"),
    get_col(base_map, "ea"),
    get_col(base_map, "oa")
  )

  base <- read_gwas(base_fp, base_cols)
  setnames(base, base_cols, c("chr", "bp", "ea", "oa"))
  base[, chr := as.character(chr)]
  base[, bp  := as.integer(bp)]
  base[, ea  := as.character(ea)]
  base[, oa  := as.character(oa)]
  base <- unique(base)

  gwas <- gwas[!base, on = c("chr", "bp", "ea", "oa")]
  cat(nrow(gwas), "novel variants after antijoin\n")
}

# add build tag
gwas[, build := input_gwas_build]

cat("Writing:", out, "\n")
data.table::fwrite(gwas, out, sep = "\t")