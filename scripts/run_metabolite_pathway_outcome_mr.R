#### DESCRIPTION #######################################
# This script runs the MR analysis of BMI on metabolite levels

#### SET INPUT #########################################
out_list     <- snakemake@input[["out_list"]]
out_map      <- snakemake@input[["out_map"]]
exp_gwas     <- snakemake@input[["exp_gwas"]]
out_maps     <- snakemake@params[["out_maps"]]
pfile        <- snakemake@params[["pfile"]]
plink2       <- NULL # on the docker path
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
cat("Reading metabolite GWAS file:", basename(exp_gwas), "\n")
m_gwas <- fread(exp_gwas)
m_gwas[, chr := as.character(chr)]
instrument_id <- m_gwas[1, id]

# create the metabolite GWAS
cat("Creating metabolite GWAS object\n")
m_g <- genepi.utils::GWAS(dat            = m_gwas,
                          map            = c(rsid = "rsid", chr = "chr", bp = "bp", ea = "ea", oa = "oa", eaf = "eaf", beta = "beta", se = "se", p = "p"),
                          fill           = TRUE,
                          drop           = TRUE,
                          fill_rsid      = FALSE,
                          missing_rsid   = "fill_CHR:BP",
                          filters        = list(),
                          id             = instrument_id,
                          trait          = instrument_id)


# cycle the outcomes files
results <- lapply(seq_along(out_list), function(i) {

  cat("Processing outcome phenotype", basename(out_list[[i]]), "\n")

  # read the gwas
  cat(" - reading GWAS\n")
  outcome_gwas <- fread(out_list[[i]])
  outcome_map <- out_maps[[i]][["map"]]
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

    # run the MR
    cat(" - running MR\n")
    res <- run_mr(mr, corr = FALSE, methods = c('mr_ivw','mr_egger','mr_weighted_median','mr_weighted_mode'))

    # add ontology mapping of metabolite
    cat(" - cleaning MR result\n")
    res <- res[, list(exposure, n_snp, fstat, outcome, method, b, b_se, p, intercept, int_se, int_p, qstat, qstat_p)]

    # log the metabolite platform
    res[, error := NA_character_]

    # return
    res

  }, error = function(e) {
    msg <- e$message
    message(msg)
    return(data.table(exposure         = instrument_id,
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