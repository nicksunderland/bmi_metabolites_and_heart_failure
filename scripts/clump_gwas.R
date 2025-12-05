#### DESCRIPTION #######################################
# This script clumps a GWAS

#### SET INPUT #########################################
gwas_file  <- snakemake@input[["gwas_file"]]
map        <- snakemake@params[["map"]]
clump_r2   <- snakemake@params[["r2"]]
clump_p1   <- snakemake@params[["p1"]]
clump_kb   <- snakemake@params[["kb"]]
pfile      <- snakemake@params[["pfile"]]
trait      <- snakemake@params[["trait"]]
id         <- snakemake@params[["id"]]
clump_file <- snakemake@output[["clump_file"]]
########################################################

# process the GWAS file
cat("Processing GWAS file:", basename(gwas_file), "\n")

# map the Yengo data to genepi.utils format
# create column map object
m <- genepi.utils::ColumnMap(lapply(map, function(x) do.call(genepi.utils::Column, x))) # nolint

# read in data to GWAS object
g <- genepi.utils::GWAS(dat            = gwas_file,
                        map            = m,
                        drop           = TRUE,
                        fill           = TRUE,
                        fill_rsid      = FALSE,
                        missing_rsid   = "fill_CHR:BP",
                        parallel_cores = parallel::detectCores(),
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

# read reference alleles
ref <- data.table::fread(paste0(pfile, ".pvar"))
ref[, "#CHROM" := as.character(`#CHROM`)]

# harmonise
gwas[ref, flip := TRUE, on = c("chr" = "#CHROM", "bp" = "POS", "oa" = "ALT", "ea" = "REF")]
gwas[flip == TRUE, c("beta", "eaf", "ea", "oa") := list(-1 * beta, 1 - eaf, oa, ea)]
gwas[!grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp, "_", oa, "_", ea)]
gwas[, flip := NULL]

# clump
clumps <- genepi.utils::clump(gwas      = gwas,
                              p1        = clump_p1,
                              r2        = clump_r2,
                              kb        = clump_kb,
                              plink2    = NULL, # should be on the PATH in the docker container
                              plink_ref = pfile)[index == TRUE, ]

# clean up
trait_param <- trait
id_param    <- id
clumps[, c("index", "clump", "trait", "id") := list(NULL, NULL, trait_param, id_param)]

# write out
data.table::fwrite(clumps, clump_file, sep = "\t")