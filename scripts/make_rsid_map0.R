#### DESCRIPTION #######################################
# This script downloads some of the metabolite GWAS files
# and makes an rsid map, to avoid doing this every time in
# the MR script

#### SET INPUT #########################################
base_file  <- snakemake@input[["base_file"]]
metab_file <- snakemake@input[["metab_file"]]
out        <- snakemake@output[["out"]]
########################################################

library(data.table)
library(fst)

# set up
setDTthreads(0)
current_threads <- getDTthreads()
cat("Current number of threads used by data.table:", current_threads, "\n")

# cols for mapping
cols <- c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele")

# get metab file
cat("Reading metabolite file\n")
metab <- read_fst(metab_file, columns = cols, as.data.table = TRUE)


if (base_file != metab_file) {
  cat("Antijoin with base file\n")
  base <- read_fst(base_file, columns = cols, as.data.table = TRUE) |> unique()
  metab <- metab[!base, on=c("chromosome", "base_pair_location", "effect_allele", "other_allele")]
  cat(nrow(metab), "additional variants\n")
}

cat("Writing variants file\n")
data.table::fwrite(metab, out, sep="\t")


