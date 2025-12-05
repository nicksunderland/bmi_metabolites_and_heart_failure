#### DESCRIPTION #######################################
# This script ...

#### SET INPUT #########################################
gwas_file  <- snakemake@input[["gwas_file"]]
gwas_build <- snakemake@params[["gwas_build"]]
gene_chr   <- snakemake@params[["gene_chr"]]
gene_start <- snakemake@params[["gene_start"]]
gene_end   <- snakemake@params[["gene_end"]]
gene_win   <- snakemake@params[["gene_win"]]
clump_r2   <- snakemake@params[["r2"]]
clump_p1   <- snakemake@params[["p1"]]
clump_kb   <- snakemake@params[["kb"]]
map        <- snakemake@params[["map"]]
fill_rsid  <- snakemake@params[["fill_rsid"]]
cores      <- snakemake@params[["cores"]]
dbsnp_dir  <- snakemake@params[["dbsnp_dir"]]
pfile      <- snakemake@params[["pfile"]]
trait      <- snakemake@params[["trait"]]
id         <- snakemake@params[["id"]]
out_file   <- snakemake@output[["out_file"]]
########################################################

library(S7)

# process the GWAS file
cat("Processing GWAS file:", basename(gwas_file), "\n")

# create column map object
cat("Creating column mapping\n")
m <- genepi.utils::ColumnMap(lapply(map, function(x) do.call(genepi.utils::Column, x))) # nolint

# subset GWAS to region
cat("Subsetting raw GWAS for gene region\n")
if (grepl("\\.(gz|tsv|csv|txt)$", gwas_file)) {
  raw_gwas    <- data.table::fread(gwas_file)
} else if (grepl("\\.(fst)$", gwas_file)) {
  raw_gwas    <- fst::read_fst(gwas_file, as.data.table = TRUE)
}
raw_chr_col <- m@map[["chr"]]@alias[1]
raw_bp_col  <- m@map[["bp"]]@alias[1]
raw_ea_col  <- m@map[["ea"]]@alias[1]
raw_oa_col  <- m@map[["oa"]]@alias[1]

# lift if needed
gwas_build <- match.arg(gwas_build, choices = c("Hg38", "Hg19"))
if (gwas_build == "Hg38") {
  raw_gwas <- genepi.utils::lift(raw_gwas, from = "Hg38", to = "Hg19", chr_col = raw_chr_col, pos_col = raw_bp_col, ea_col = raw_ea_col, oa_col = raw_oa_col, remove_duplicates = FALSE)
}

# subset on build 37 coordinates
raw_gwas    <- raw_gwas[get(raw_chr_col) == gene_chr & get(raw_bp_col) > gene_start - gene_win & get(raw_bp_col) < gene_end + gene_win]

# read in data to GWAS object
cat("Creating GWAS object for gene region\n")
g <- genepi.utils::GWAS(dat            = raw_gwas,
                        map            = m,
                        drop           = TRUE,
                        fill           = TRUE,
                        fill_rsid      = fill_rsid,
                        missing_rsid   = "fill_CHR:BP",
                        parallel_cores = cores,
                        dbsnp_dir      = dbsnp_dir,
                        trait          = trait,
                        id             = id,
                        filters        = list(
                          beta_invalid    = "!is.infinite(beta) & abs(beta) < 20",
                          eaf_invalid     = "eaf > 0 & eaf < 1",
                          p_invalid       = "!is.infinite(p)",
                          se_invalid      = "!is.infinite(se)",
                          alleles_invalid = "!is.na(ea) & !is.na(oa)",
                          chr_missing     = "!is.na(chr)",
                          bp_missing      = "!is.na(bp)",
                          beta_missing    = "!is.na(beta)",
                          se_missing      = "!is.na(se)",
                          p_missing       = "!is.na(p)",
                          eaf_missing     = "!is.na(eaf)"),
                        verbose      = TRUE)

# to data.table
gwas <- genepi.utils::as.data.table(g)

# # read reference alleles
# ref <- data.table::fread(paste0(pfile, ".pvar"))
# ref[, "#CHROM" := as.character(`#CHROM`)]
#
# # harmonise
# gwas[ref, flip := TRUE, on = c("chr" = "#CHROM", "bp" = "POS", "oa" = "ALT", "ea" = "REF")]
# gwas[flip == TRUE, c("beta", "eaf", "ea", "oa") := list(-1 * beta, 1 - eaf, oa, ea)]
# gwas[!grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp, "_", oa, "_", ea)]
# gwas[, flip := NULL]

# clump
cat("Clumping GWAS object\n")
clumps <- genepi.utils::clump(gwas      = gwas,
                              p1        = clump_p1,
                              r2        = clump_r2,
                              kb        = clump_kb,
                              plink2    = NULL, # should be on the PATH in the docker container
                              plink_ref = pfile)[index == TRUE, ]

# clean up
clumps[, c("index", "clump") := NULL]

# error if no instrument
if (nrow(clumps) == 0) {
  cat("Error - no variants found, check the gene region and clumping parameters\n")
  stop("Invalid instrument produced")
}

# write out
cat("Saving file\n")
data.table::fwrite(clumps, out_file, sep = "\t")