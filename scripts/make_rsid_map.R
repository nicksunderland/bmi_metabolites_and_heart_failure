#### DESCRIPTION #######################################
# This script downloads some of the metabolite GWAS files
# and makes an rsid map, to avoid doing this every time in
# the MR script

#### SET INPUT #########################################
metab_files <- snakemake@input[["metab_files"]]
dbsnp_dir   <- snakemake@params[["dbsnp_dir"]]
cores       <- snakemake@resources[["cpus_per_task"]]
platform    <- snakemake@params[["platform"]]
map_file    <- snakemake@output[["map_file"]]
########################################################

library(data.table)
library(fst)

# set up
setDTthreads(0)
current_threads <- getDTthreads()
cat("Current number of threads used by data.table:", current_threads, "\n")
cat("There are", length(metab_files), "metabolite gwas files to read\n")

# cols for mapping
map <- NULL
for (i in seq_along(metab_files)) {
  f <- metab_files[[i]]
  cat("Processing:", basename(f), "(", i, "/", length(metab_files), ")\n")
  d <- fread(f)

  if (is.null(map)) {
    map <- d
  } else {
    map <- rbind(map, d)
  }

  # Apply unique every 20 files or on the last iteration
  if (i %% 50 == 0 || i == length(metab_files)) {
    map <- unique(map)
  }
}

map <- unique(map)


# annotate rsid
map <- genepi.utils::chrpos_to_rsid(map,
                      chr_col = "chromosome",
                      pos_col = "base_pair_location",
                      ea_col  = "effect_allele",
                      nea_col = "other_allele",
                      flip    = "allow",
                      build   = ifelse(platform == "nightingale","b37_dbsnp156","b38_dbsnp156"),
                      parallel_cores = cores,
                      dbsnp_dir = dbsnp_dir)

data.table::setnames(map, "RSID", ifelse(platform == "nightingale","rsid_b37","rsid_b38"))
pos_col <- ifelse(platform == "nightingale","base_pair_location_b37","base_pair_location_b38")
map[, (pos_col) := base_pair_location]

# liftover
map <- genepi.utils::lift(map,
                          from = ifelse(platform == "nightingale","Hg19","Hg38"),
                          to   = ifelse(platform == "nightingale","Hg38","Hg19"),
                          chr_col = "chromosome",
                          pos_col = "base_pair_location",
                          ea_col  = "effect_allele",
                          oa_col  = "other_allele",
                          remove_duplicates = FALSE)

# annotate rsid other way
map <- genepi.utils::chrpos_to_rsid(map,
                      chr_col = "chromosome",
                      pos_col = "base_pair_location",
                      ea_col  = "effect_allele",
                      nea_col = "other_allele",
                      flip    = "allow",
                      build   = ifelse(platform == "nightingale","b38_dbsnp156","b37_dbsnp156"),
                      parallel_cores = cores,
                      dbsnp_dir = dbsnp_dir)

data.table::setnames(map, "RSID", ifelse(platform == "nightingale","rsid_b38","rsid_b37"))
data.table::setnames(map, "base_pair_location", ifelse(platform == "nightingale","base_pair_location_b38","base_pair_location_b37"))


write_fst(map, map_file, compress=100)


