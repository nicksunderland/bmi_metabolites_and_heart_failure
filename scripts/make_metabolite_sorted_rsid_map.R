#### DESCRIPTION #######################################
# This script creates a mapping file from the original
# metabolite mapping file, adding a sorted ID
#
#### SET INPUT #########################################
orig_map  <- snakemake@input[["orig_map"]]
dbsnp_dir <- snakemake@params[["dbsnp_dir"]]
cores     <- snakemake@resources[["cpus_per_task"]]
map_out   <- snakemake@output[["map_out"]]
########################################################

# requirements
library(data.table)
library(fst)

# set up
setDTthreads(0)
current_threads <- getDTthreads()
cat("Current number of threads used by data.table:", current_threads, "\n")


# cols for mapping
cat("Reading original mapping file\n")
map <- read_fst(orig_map, as.data.table = TRUE)
map[, variant_id := NULL]
setnames(map, c("chromosome","base_pair_location_b37","base_pair_location_b38","effect_allele","other_allele"), c("chr","bp_b37","bp_b38","ea", "oa"))


# get the ok rsids
cat("Assessing b37 RSID coding\n")
ok  <- map[grepl("^rs[0-9]+$", rsid_b37)]
bad <- map[!grepl("^rs[0-9]+$", rsid_b37)]

# try re-annotate rsid
cat("Reannotating", nrow(bad), "badly formatted b37 RSIDs\n")
bad[, rsid_b37 := NULL]
bad <- genepi.utils::chrpos_to_rsid(bad,
                      chr_col = "chr",
                      pos_col = "bp_b37",
                      ea_col  = "ea",
                      nea_col = "oa",
                      flip    = "allow",
                      build   = "b37_dbsnp156",
                      parallel_cores = cores,
                      dbsnp_dir = dbsnp_dir)
data.table::setnames(bad, "RSID","rsid_b37")

# set missing rsid to chr:pos
cat("Setting", bad[!grepl("^rs[0-9]+$", rsid_b37), .N], "remaining badly formatted b37 RSIDs to chr:bp\n")
bad[!grepl("^rs[0-9]+$", rsid_b37), rsid_b37 := paste0(chr, ":", bp_b37)]

# rbind
map <- rbind(ok, bad)


# get the ok rsids
cat("Assessing b38 RSID coding\n")
ok  <- map[grepl("^rs[0-9]+$", rsid_b38)]
bad <- map[!grepl("^rs[0-9]+$", rsid_b38)]

# try re-annotate rsid
cat("Reannotating", nrow(bad), "badly formatted b38 RSIDs\n")
bad[, rsid_b38 := NULL]
bad <- genepi.utils::chrpos_to_rsid(bad,
                      chr_col = "chr",
                      pos_col = "bp_b38",
                      ea_col  = "ea",
                      nea_col = "oa",
                      flip    = "allow",
                      build   = "b38_dbsnp156",
                      parallel_cores = cores,
                      dbsnp_dir = dbsnp_dir)
data.table::setnames(bad, "RSID","rsid_b38")

# set missing rsid to chr:pos
cat("Setting", bad[!grepl("^rs[0-9]+$", rsid_b38), .N], "remaining badly formatted b38 RSIDs to chr:bp\n")
bad[!grepl("^rs[0-9]+$", rsid_b38), rsid_b38 := paste0(chr, ":", bp_b38)]

# rbind
map <- rbind(ok, bad)


# create sorted allele id
cat("Creating sorted allele variant ID\n")
map[, variant_id_b37 := ifelse(ea < oa, paste(rsid_b37, ea, oa, sep="_"), paste(rsid_b37, oa, ea, sep="_"))]
map[, variant_id_b38 := ifelse(ea < oa, paste(rsid_b38, ea, oa, sep="_"), paste(rsid_b38, oa, ea, sep="_"))]
setcolorder(map, c("variant_id_b37", "variant_id_b38"))

# save
cat("Saving to .fst\n")
map <- unique(map)
setkeyv(map, c("chr", "bp_b37"))
write_fst(map, map_out, compress=100)


