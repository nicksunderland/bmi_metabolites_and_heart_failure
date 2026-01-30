#### DESCRIPTION #######################################
# This script creates a mapping file from the hermes
# incidence GWASs and adds a sorted ID
#
#### SET INPUT #########################################
gwases    <- snakemake@input[["gwases"]]
dbsnp_dir <- snakemake@params[["dbsnp_dir"]]
gwas_maps <- snakemake@params[["gwas_maps"]]
cores     <- snakemake@resources[["cpus_per_task"]]
map_out   <- snakemake@output[["map_out"]]
########################################################

library(data.table)
library(fst)

# set up
setDTthreads(0)
current_threads <- getDTthreads()
cat("Current number of threads used by data.table:", current_threads, "\n")
cat("There are", length(gwases), "outcome gwas files to read\n")
print(str(gwases))
print(str(gwas_maps))

# cols for mapping
map <- NULL
for (i in seq_along(gwases)) {
  f <- gwases[[i]]
  cat("Processing file:", basename(f), "(", i, "/", length(gwases), ")\n")
  gwas_map <- gwas_maps[[i]]
  old_names <- sapply(c("rsid", "chr", "bp", "ea", "oa"), function(x) gwas_map[[x]][["alias"]])
  d <- fread(f, select = unlist(unname(old_names)), col.names = names(old_names))
  d[, `:=`(chr = as.character(chr), ea = toupper(ea), oa = toupper(oa))]
  if (is.null(map)) {
    map <- d
  } else {
    map <- unique(rbind(map, d))
  }
}

# get the ok rsids
cat("Assessing RSID coding\n")
ok  <- map[grepl("^rs[0-9]+$", rsid)]
bad <- map[!grepl("^rs[0-9]+$", rsid)]

if (nrow(bad) > 0) {
  # try re-annotate rsid
  cat("Reannotating", nrow(bad), "badly formatted RSIDs\n")
  bad[, rsid := NULL]
  bad <- genepi.utils::chrpos_to_rsid(bad,
                        chr_col = "chr",
                        pos_col = "bp",
                        ea_col  = "ea",
                        nea_col = "oa",
                        flip    = "allow",
                        build   = "b37_dbsnp156",
                        parallel_cores = cores,
                        dbsnp_dir = dbsnp_dir)
  data.table::setnames(bad, "RSID","rsid")

  # set missing rsid to chr:pos
  cat("Setting", bad[!grepl("^rs[0-9]+$", rsid), .N], "remaining badly formatted RSIDs to chr:bp\n")
  bad[!grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]

  # rbind
  map <- rbind(ok, bad[, .(rsid, chr, bp, ea, oa)])

} else {
  map <- ok
}

# create sorted allele id
cat("Creating sorted allele variant ID\n")
map[, variant_id := ifelse(ea < oa, paste(rsid, ea, oa, sep="_"), paste(rsid, oa, ea, sep="_"))]
setcolorder(map, "variant_id")

# save
cat("Saving to .fst\n")
map[, chr := as.character(chr)]
map <- map[!is.na(chr) & chr != ""]
map <- map[!is.na(bp)]
map <- map[!is.na(ea)]
map <- map[!is.na(oa)]
map <- unique(map)
setkeyv(map, c("chr", "bp"))
write_fst(map, map_out, compress=100)


