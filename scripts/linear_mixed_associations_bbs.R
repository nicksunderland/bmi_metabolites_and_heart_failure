#### DESCRIPTION #######################################
# Purpose of script:
# Take in the processed metabolite data.table and run
# linear mixed model association analyses.
#
# Author: Nick Sunderland
#
# Date Created: 2025-04-08
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
metaoprep_obj <- snakemake@input[["metaoprep_obj"]]
model_df      <- snakemake@input[["model_df"]]
functions     <- snakemake@params[["functions"]]
data_type     <- snakemake@params[["data_type"]]
platform      <- snakemake@params[["platform"]]
assoc         <- snakemake@output[["assoc"]]
log_file      <- snakemake@output[["log_file"]]
########################################################


# requirements
suppressPackageStartupMessages(library(data.table))
library(cli)
library(future)
library(furrr)
library(progressr)
source(functions)

# logging
sink(log_file); sink()


# start
log("By-Band-Sleeve Linear Mixed Model Association Analysis", func = "cli_h1")
log("Date:", as.character(Sys.Date()), func = "cli_text")


# read object
log("Reading in model data...", func = "cli_alert_info")
model_dat <- fread(model_df)


# clean factors
log("Cleaning columns factors...", func = "cli_alert_info")
model_dat[, `:=`(timepoint = factor(timepoint, levels = c("baseline", "end")),
                 sex       = factor(sex, labels = c("male", "female")),
                 study     = factor(study, levels = c("bbs_mainstudy", "bbs_substudy")),
                 site      = as.factor(site),
                 study_id  = as.factor(study_id))]


# run linear mixed model
log("Running linear mixed model...", func = "cli_alert_info")
metabs <- grep("^metab_", names(model_dat), value=TRUE)
plan(multisession, workers = parallel::detectCores()-1)
with_progress({
  p    <- progressor(steps = length(metabs))
  res  <- future_map(stats::setNames(metabs, metabs), run_model,
                     model_dat = model_dat,
                     ns_func   = "lmerTest::lmer",
                     fixed     = c("timepoint", "age", "sex", "site", "study", "time_in_freezer"),
                     random    = "study_id",
                     main_study= model_dat$main_study[1],
                     data_type = data_type,
                     model_name= "linear_mixed_time",
                     platform  = platform,
                     p         = p)

})
res <- rbindlist(res, idcol = "feature_id", fill=TRUE)
res[, p.value_holm := stats::p.adjust(p.value, method="holm", n=.N), by="term"]
res[, p.value_fdr := stats::p.adjust(p.value, method="fdr", n=.N), by="term"]
log("Linear mixed models complete.", func = "cli_alert_success")


# annotate metabolite names
log("Annotating metabolites...", func = "cli_alert_info")
m          <- readRDS(metaoprep_obj)
features   <- m$features
annot_cols <- c("feature_id", "label", "raw_id", "pathway_group", "pathway", "sub_pathway", "derived_feature", "missingness", "outlier_count", "independent", "independent_k")
res        <- res[features[, .SD, .SDcols = annot_cols], on="feature_id"]
setcolorder(res, annot_cols)


# write out
log("Writing results...", func = "cli_alert_info")
fwrite(res, assoc, sep = "\t")


# capture session info
log("\n", func = "cli_text")
log("Session information:", func = "cli_h2")
si <- capture.output(sessionInfo())
for (i in si) log(i, func = "cli_text")
