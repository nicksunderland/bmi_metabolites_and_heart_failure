#### DESCRIPTION #######################################
# Merges per-GWAS variant extracts, annotates RSIDs in
# both builds via liftover, and writes the global variant map

#### SET INPUT #########################################
variant_files <- snakemake@input[["variant_files"]]
dbsnp_dir     <- snakemake@params[["dbsnp_dir"]]
cores         <- as.integer(snakemake@resources[["cpus_per_task"]])
map_file      <- snakemake@output[["map_file"]]
########################################################

library(data.table)
library(fst)

setDTthreads(cores)
cat("Threads:", getDTthreads(), "\n")
cat("Merging", length(variant_files), "variant files\n")

# stack variant files into separate build tables
hg19_map <- NULL
hg38_map <- NULL

for (i in seq_along(variant_files)) {
  f <- variant_files[[i]]
  cat("Reading:", basename(f), "(", i, "/", length(variant_files), ")\n")
  d <- fread(f)
  d[, ea := toupper(ea)]
  d[, oa := toupper(oa)]
  build <- d[1, build]
  if (build == "Hg19") {
    hg19_map <- if (is.null(hg19_map)) d else rbind(hg19_map, d, fill = TRUE)
    if (i %% 50 == 0) hg19_map <- unique(hg19_map, by = c("chr", "bp", "ea", "oa"))
  } else {
    hg38_map <- if (is.null(hg38_map)) d else rbind(hg38_map, d, fill = TRUE)
    if (i %% 50 == 0) hg38_map <- unique(hg38_map, by = c("chr", "bp", "ea", "oa"))
  }
}

hg19_map <- if (!is.null(hg19_map)) unique(hg19_map, by = c("chr", "bp", "ea", "oa")) else data.table()
hg38_map <- if (!is.null(hg38_map)) unique(hg38_map, by = c("chr", "bp", "ea", "oa")) else data.table()
cat(nrow(hg19_map), "unique Hg19 variants\n")
cat(nrow(hg38_map), "unique Hg38 variants\n")


# liftover Hg38 to Hg19 (b37) only - no annotation yet
if (nrow(hg38_map) > 0) {
  cat("Lifting Hg38 variants to Hg19\n")
  hg38_map[, bp_b38 := bp]
  hg38_map <- genepi.utils::lift(hg38_map,
    from    = "Hg38",
    to      = "Hg19",
    chr_col = "chr",
    pos_col = "bp",
    ea_col  = "ea",
    oa_col  = "oa",
    remove_duplicates = TRUE)
}

if (nrow(hg19_map) > 0) {
  hg19_map[, bp_b38 := NA_integer_]
}

# combine
map <- rbind(hg19_map, hg38_map, fill = TRUE)
map[, build := NULL]
map <- unique(map, by = c("chr", "bp", "bp_b38", "ea", "oa"))
cat(nrow(map), "unique variants on b37 coordinates\n")

map <- genepi.utils::chrpos_to_rsid(map[!is.na(bp)],
  chr_col        = "chr",
  pos_col        = "bp",
  ea_col         = "ea",
  nea_col        = "oa",
  flip           = "allow",
  build          = "b37_dbsnp156",
  parallel_cores = cores,
  dbsnp_dir      = dbsnp_dir)
setnames(map, "RSID", "rsid")

# fallback rsid formatting
cat("Setting", map[!grepl("^rs[0-9]+$", rsid), .N], "bad RSIDs to chr:bp\n")
map[!grepl("^rs[0-9]+$", rsid), rsid := paste0(chr, ":", bp)]

# variant IDs
map[, variant_id := ifelse(ea < oa,
  paste(rsid, ea, oa, sep = "_"),
  paste(rsid, oa, ea, sep = "_"))]

setcolorder(map, c("variant_id", "rsid", "chr", "bp", "ea", "oa", "bp_b38"))

cat("Writing map:", map_file, "\n")
fwrite(map, map_file, sep="\t")