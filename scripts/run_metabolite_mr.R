#### DESCRIPTION #######################################
# This script runs the MR analysis of BMI on metabolite levels

#### SET INPUT #########################################
exp_path       <- snakemake@input[["exp_path"]]
rsid_map       <- snakemake@input[["rsid_map"]]
met_path       <- snakemake@input[["met_path"]]
met_yaml       <- snakemake@input[["met_yaml"]]
platform       <- snakemake@params[["platform"]]
exposure       <- snakemake@params[["exposure"]]
metabolite_id  <- snakemake@params[["metabolite_id"]]
exp_clump_map  <- snakemake@params[["exp_clump_map"]]
metabolite_map <- snakemake@params[["metabolite_map"]]
pfile          <- snakemake@params[["pfile"]]
plink2         <- NULL
results_file   <- snakemake@output[["results_file"]]
########################################################

library(S7)
library(data.table)
library(fst)

cat("Running MR analyses:", exposure, "vs.", metabolite_id, "\n")


# create exposure map
cat("[i] Creating", exposure, "file mapping\n")
exp_map <- genepi.utils::ColumnMap(lapply(exp_clump_map, function(x) do.call(genepi.utils::Column, x)))


# get unique chromosome in exposure (to use to speed up processing of metab gwas in case of gene region, below)
raw_exp       <- data.table::fread(exp_path)
raw_chr_col   <- exp_map@map[["chr"]]@alias[1]
unique_chroms <- unique(raw_exp[, get(raw_chr_col)])


# create the exposure GWAS object
cat("[i] Creating", exposure, "GWAS object\n")
exp <- genepi.utils::GWAS(dat          = raw_exp,
                          map          = exp_map,
                          fill         = TRUE,
                          drop         = TRUE,
                          fill_rsid    = FALSE,
                          missing_rsid = "fill_CHR:BP",
                          id           = exposure,
                          trait        = exposure)


# read metabolite GWAS
cat("[i] Reading metabolite GWAS data\n")
metab_dat <- read_fst(met_path, as.data.table = TRUE)
metab_dat <- unique(metab_dat, by=c("chromosome", "base_pair_location", "effect_allele", "other_allele"))
meta      <- yaml::read_yaml(met_yaml)


# subset for chromosomes in exposure to speed up GWAS/MR object creation
metab_dat[, chromosome := as.character(chromosome)]
metab_dat <- metab_dat[chromosome %in% as.character(unique_chroms)]


# add the rsids and map to b37
cat("[i] Merging with RSID map\n")
rs_fst <- fst(rsid_map)
rs_chr <- as.character(rs_fst[, "chromosome"])

if (platform == "nightingale") {

  rs_map <- as.data.table(rs_fst[which(rs_chr %in% as.character(unique_chroms)), c("rsid_b37", "chromosome", "base_pair_location_b37", "effect_allele", "other_allele")])
  setnames(rs_map, c("rsid_b37", "base_pair_location_b37"), c("rsid", "base_pair_location"))
  rs_map <- unique(rs_map, by=c("chromosome", "base_pair_location", "effect_allele", "other_allele"))
  metab_dat[rs_map, rsid := i.rsid, on = c("chromosome", "base_pair_location", "effect_allele", "other_allele")]

} else if (platform == "metabolon") {

  rs_map <- as.data.table(rs_fst[which(rs_chr %in% as.character(unique_chroms)), c("rsid_b37", "chromosome", "base_pair_location_b37", "base_pair_location_b38", "effect_allele", "other_allele")])
  rs_map <- unique(rs_map, by=c("chromosome", "base_pair_location_b38", "effect_allele", "other_allele"))
  metab_dat[rs_map, c("rsid", "base_pair_location") := list(i.rsid_b37, i.base_pair_location_b37),
            on = c("chromosome"="chromosome", "base_pair_location"="base_pair_location_b38", "effect_allele"="effect_allele", "other_allele"="other_allele")]
  metab_dat <- metab_dat[!is.na(base_pair_location)]

}


# map for metabolite GWAS
cat("[i] Creating metabolite file mapping\n")
metabolite_map[["rsid"]] <- list(name = "rsid", alias = "rsid", type =  "character") # update column
metab_map <- genepi.utils::ColumnMap(lapply(metabolite_map, function(x) do.call(genepi.utils::Column, x)))


# create the metabolite GWAS object
cat("[i] Creating metabolite GWAS object\n")
metab <- genepi.utils::GWAS(dat            = metab_dat,
                            map            = metab_map,
                            fill           = TRUE,
                            drop           = TRUE,
                            fill_rsid      = FALSE,
                            missing_rsid   = "fill_CHR:BP",
                            id             = metabolite_id,
                            trait          = meta$trait_description,
                            n              = sum(sapply(meta$samples, function(x) x[["sample_size"]])))


# get proxies and subset
cat(" - looking for proxies for missing variants in outcome\n")
metab <- genepi.utils::get_proxies(metab,
                                   snps       = exp@rsid,
                                   then       = "subset",
                                   stat       = "r2-unphased",
                                   win_kb     = 125,
                                   win_r2     = 0.8,
                                   win_ninter = Inf,
                                   proxy_eaf  = NULL,
                                   plink2     = plink2,
                                   pfile      = pfile)


# create an MR object (harmonises exp and metab)
cat("[i] Creating MR object\n")
mr <- genepi.utils::MR(exp, metab, harmonise_strictness = 2)

# run the MR
cat("[i] Running MR\n")
res <- genepi.utils::run_mr(mr, corr = FALSE, methods = c('mr_ivw','mr_egger','mr_weighted_median','mr_weighted_mode'))

# add ontology mapping of metabolite
cat("[i] Cleaning MR result\n")
res <- res[, list(exposure, n_snp, fstat, outcome, method, b, b_se, p, intercept, int_se, int_p, qstat, qstat_p)]

# log the metabolite platform
res[, metabolite_platform := platform]
res[, metabolite_id := metabolite_id]
setcolorder(res, c("metabolite_platform", "metabolite_id"), after = "outcome")

# save the result
data.table::fwrite(res, results_file, sep="\t")