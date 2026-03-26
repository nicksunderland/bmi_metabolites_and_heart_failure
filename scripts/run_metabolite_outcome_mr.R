#### DESCRIPTION #######################################
# This script runs the MR analysis of BMI on metabolite levels

#### SET INPUT #########################################
out_list     <- snakemake@input[["out_list"]]
out_map      <- snakemake@input[["out_map"]]
metab_gwas   <- snakemake@input[["metab_gwas"]]
metab_meta   <- snakemake@input[["metab_meta"]]
metab_map    <- snakemake@input[["metab_map"]]
cores        <- snakemake@resources[["cpus_per_task"]]
out_maps     <- snakemake@params[["out_maps"]]
platform     <- snakemake@params[["platform"]]
clump_r2     <- snakemake@params[["clump_r2"]]
clump_p1     <- snakemake@params[["clump_p1"]]
clump_kb     <- snakemake@params[["clump_kb"]]
pfile        <- snakemake@params[["pfile"]]
plink2       <- NULL # on the docker path
mr_instrument<- snakemake@output[["mr_instrument"]]
mr_result    <- snakemake@output[["mr_result"]]
########################################################

# requirements
library(S7)
library(genepi.utils)
library(data.table)
library(fst)
library(yaml)
library(parallel)


# set up
setDTthreads(0)
current_threads <- getDTthreads()
cat("Current number of threads used by data.table:", current_threads, " \n")
cat("There are", length(out_list), "OUTCOME gwas files to process\n")


# read the HERMES map
cat("Reading OUTCOME map\n")
o_map <- read_fst(out_map, as.data.table = TRUE)


# read the metabolite gwas
cat("Reading metabolite GWAS file:", basename(metab_gwas), "\n")
m_gwas <- read_fst(metab_gwas, as.data.table = TRUE)
m_gwas[, chromosome := as.character(chromosome)]

# read the metabolite meta data
cat("Reading metabolite meta-data file:", basename(metab_meta), "\n")
m_meta <- read_yaml(metab_meta)


# read the metabolite meta data
cat("Reading metabolite map file:", basename(metab_map), "\n")
m_map <- read_fst(metab_map, as.data.table = TRUE)


# add the variant ID from mapping
cat("Annotating metabolite GWAS with mapping variant id\n")
if (platform == "nightingale") {
  m_gwas[m_map, `:=`(variant_id = i.variant_id_b37, rsid = i.rsid_b37), on = c("chromosome"="chr", "base_pair_location"="bp_b37", "effect_allele"="ea", "other_allele"="oa")]
} else if (platform == "metabolon") {
  m_gwas[m_map, `:=`(variant_id = i.variant_id_b37, rsid = i.rsid_b37, base_pair_location = i.bp_b37), on = c("chromosome"="chr", "base_pair_location"="bp_b38", "effect_allele"="ea", "other_allele"="oa")]
}


# clump metabolite GWAS
setnames(m_gwas, c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"), c("rsid", "chr", "bp", "ea", "oa", "eaf", "beta", "se", "p"))
clumps <- genepi.utils::clump(gwas      = m_gwas,
                              p1        = clump_p1,
                              r2        = clump_r2,
                              kb        = clump_kb,
                              parallel_cores = cores,
                              plink2    = plink2, # should be on the PATH in the docker container
                              plink_ref = pfile)[index == TRUE, ]


# create the metabolite GWAS
cat("Creating metabolite GWAS object\n")
m_g <- genepi.utils::GWAS(dat            = clumps,
                          map            = c(rsid = "rsid", chr = "chr", bp = "bp", ea = "ea", oa = "oa", eaf = "eaf", beta = "beta", se = "se", p = "p"),
                          fill           = TRUE,
                          drop           = TRUE,
                          fill_rsid      = FALSE,
                          missing_rsid   = "fill_CHR:BP",
                          filters        = list(),
                          id             = m_meta$gwas_id,
                          trait          = m_meta$trait_description,
                          n              = sum(sapply(m_meta$samples, function(x) x[["sample_size"]])))


# metabolite instrument
cat("Storing metabolite instrument\n")
instrument <- genepi.utils::as.data.table(m_g)
print(head(instrument, 5))


# cycle the outcomes files
results <- lapply(seq_along(out_list), function(i) {

  cat("Processing outcome phenotype", basename(out_list[[i]]), "\n")

  # update instrument df
  used_col  <- paste0("included_", names(out_maps)[i])
  proxy_col <- paste0("proxy_", names(out_maps)[i])
  beta_col  <- paste0("beta_", names(out_maps)[i])
  se_col    <- paste0("se_", names(out_maps)[i])
  p_col     <- paste0("p_", names(out_maps)[i])
  instrument[, c(used_col, proxy_col, beta_col, se_col, p_col) := list(FALSE, NA_character_, NA_real_, NA_real_, NA_real_)]

  # read the gwas
  cat(" - reading GWAS\n")
  outcome_gwas <- fread(out_list[[i]])
  outcome_map <- out_maps[[i]]
  old_names <- sapply(outcome_map, function(x) x[["alias"]])
  new_names <- sapply(outcome_map, function(x) x[["name"]])
  setnames(outcome_gwas, old_names, new_names)
  outcome_gwas[, `:=`(rsid = as.character(rsid), chr = as.character(chr), ea = toupper(ea), oa = toupper(oa))]

  # join with the map to update rsids for HERMES
  cat(" - joining with map RSIDs\n")
  outcome_gwas[o_map, `:=`(rsid = i.rsid), on = c("chr", "bp", "ea", "oa")]

  res <- tryCatch({

    # create the GWAS object
    cat(" - creating OUTSOME GWAS object\n")
    o_g <- genepi.utils::GWAS(dat            = outcome_gwas,
                              map            = c(rsid = "rsid", chr = "chr", bp = "bp", ea = "ea", oa = "oa", eaf = "eaf", beta = "beta", se = "se", p = "p", n = "n", ncase = "ncase"),
                              fill           = TRUE,
                              drop           = TRUE,
                              fill_rsid      = FALSE,
                              missing_rsid   = "fill_CHR:BP",
                              filters        = list(),
                              id             = basename(out_list[[i]]),
                              trait          = names(out_maps)[i])

    # get proxies and subset
    cat(" - looking for proxies for missing variants\n")
    o_g <- get_proxies(o_g,
                       snps       = m_g@rsid,
                       then       = "subset",
                       stat       = "r2-unphased",
                       win_kb     = 125,
                       win_r2     = 0.8,
                       win_ninter = Inf,
                       proxy_eaf  = NULL,
                       plink2     = plink2,
                       pfile      = pfile)

    # make MR object
    cat(" - creating MR object\n")
    if (length(m_g@rsid) < 5) {
      cat("[i] harmonising without excluding palindromic SNPs with indeterminate freq to maintain instrument size\n")
      strictness <- 1
    } else {
      strictness <- 2
    }
    mr <- MR(m_g, o_g, harmonise_strictness = strictness)

    # get the SNPs used in the instrument
    used_instrument <- genepi.utils::as.data.table(mr)
    instrument[used_instrument, c(used_col, proxy_col, beta_col, se_col, p_col) :=  list(TRUE, i.proxy_o_snp, i.by, i.byse, i.py), on = c("chr","bp","ea","oa")]
    instrument[used_instrument, c(used_col, proxy_col, beta_col, se_col, p_col) :=  list(TRUE, i.proxy_o_snp, -1*i.by, i.byse, i.py), on = c("chr"="chr","bp"="bp","ea"="oa","oa"="ea")]

    # run the MR
    cat(" - running MR\n")
    res <- run_mr(mr, corr = FALSE, methods = c('mr_ivw','mr_egger','mr_weighted_median','mr_weighted_mode'))

    # add ontology mapping of metabolite
    cat(" - cleaning MR result\n")
    res <- res[, list(exposure, n_snp, fstat, outcome, method, b, b_se, p, intercept, int_se, int_p, qstat, qstat_p)]

    # log the metabolite platform
    res[, metabolite_platform := platform]
    res[, metabolite_id := m_meta$gwas_id]
    res[, error := NA_character_]
    setcolorder(res, c("metabolite_platform", "metabolite_id"), after = "exposure")

    # return
    res

  }, error = function(e) {
    msg <- e$message
    message(msg)
    return(data.table(exposure         = m_meta$trait_description,
                      metabolite_platform = platform,
                      metabolite_id    = m_meta$gwas_id,
                      n_snp            = NA_integer_,
                      fstat            = NA_real_,
                      outcome          = names(out_maps)[i],
                      method           = NA_character_,
                      b                = NA_real_,
                      b_se             = NA_real_,
                      p                = NA_real_,
                      intercept        = NA_real_,
                      int_se           = NA_real_,
                      int_p            = NA_real_,
                      qstat            = NA_real_,
                      qstat_p          = NA_real_,
                      error            = msg))
  })

  # return
  res

}) |> rbindlist(fill=TRUE)


# save the result
data.table::fwrite(results, mr_result, sep="\t")
data.table::fwrite(instrument, mr_instrument, sep="\t")