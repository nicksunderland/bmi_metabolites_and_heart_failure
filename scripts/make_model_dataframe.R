#### DESCRIPTION #######################################
# Purpose of script:
# Take in the processed metabolite data in the form of a
# metaboprep2::Metabolites object and wrangle the data into
# a data.table suitable for running the association analyses.
#
# Author: Nick Sunderland
#
# Date Created: 2025-04-08
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
qc_data_obj <- snakemake@input[["qc_data_obj"]]
functions   <- snakemake@params[["functions"]]
platform    <- snakemake@params[["platform"]]
main_study  <- snakemake@params[["study"]]
model_raw   <- snakemake@output[["model_raw"]]
model_raw_z <- snakemake@output[["model_raw_z"]]
model_rnt   <- snakemake@output[["model_rnt"]]
log_file    <- snakemake@output[["log_file"]]
########################################################


# requirements
suppressPackageStartupMessages(library(data.table))
library(cli)


# helper functions
source(functions)


# logging
sink(log_file); sink()


# start
log("By-Band-Sleeve Linear Model Data Wrangling", func = "cli_h1")
log("Date:", as.character(Sys.Date()), func = "cli_text")


# read object
log("Reading in metaboprep2 object...", func = "cli_alert_info")
m <- readRDS(qc_data_obj)


# get data
log("Extracting required data...", func = "cli_alert_info")
m_dat    <- m$data
samples  <- m$samples
features <- m$features


# helper functions
z_score <- function(x, colname) {
  r <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  if (!startsWith(colname, "olink")) {
    r[is.na(r)] <- min(r, na.rm = TRUE)
  }
  r[is.infinite(r)] <- NA_real_
  r
}
rnt <- function(x, colname) {
  r <- rank(x, na.last = "keep", ties.method = "random")
  r <- qnorm((r - 0.5) / sum(!is.na(x)))
  if (!startsWith(colname, "olink")) {
    r[is.na(r)] <- min(r, na.rm = TRUE)
  }
  r[is.infinite(r)] <- NA_real_
  r
}


# data storage
data_list <- list(raw          = list(path = model_raw, data = NA),
                  raw_z        = list(path = model_raw_z, data = NA),
                  rnt          = list(path = model_rnt, data = NA))


# raw data
data_list[["raw"]][["data"]] <- m_dat


# raw Z score
log("Applying Z-score transformation (raw-Z)...", func = "cli_alert_info")
data_list[["raw_z"]][["data"]] <- mapply(FUN     = z_score,
                                         x       = as.data.frame(data_list[["raw"]][["data"]]),
                                         colname = colnames(data_list[["raw"]][["data"]]))
rownames(data_list[["raw_z"]][["data"]]) <- rownames(data_list[["raw"]][["data"]])


# apply rank-based inverse normal transformation
log("Applying rank inverse normal transformation (RNT)...", func = "cli_alert_info")
data_list[["rnt"]][["data"]] <- mapply(FUN     = rnt,
                                       x       = as.data.frame(data_list[["raw"]][["data"]]),
                                       colname = colnames(data_list[["raw"]][["data"]]))
rownames(data_list[["rnt"]][["data"]]) <- rownames(data_list[["raw"]][["data"]])


# merge data for model
log("Formatting data for model...", func = "cli_alert_info")
for (i in seq_along(data_list)) {
  model_dat <- data.table::as.data.table(data_list[[i]][["data"]], keep.rownames = TRUE)
  setnames(model_dat, "rn", "sample_id")
  model_dat <- samples[model_dat, on = "sample_id"]
  model_dat[, `:=`(main_study = main_study,
                   timepoint  = factor(timepoint, levels = c("Baseline", "end"), labels = c("baseline", "end")),
                   sex        = factor(sex, levels = 1:2, labels = c("male", "female")))]
  data_list[[i]][["data"]] <- model_dat
}


# write out
for (i in seq_along(data_list)) {
  log("Writing", names(data_list)[i], "data:", data_list[[i]][["path"]], func = "cli_alert_info")
  fwrite(data_list[[i]][["data"]], data_list[[i]][["path"]], sep = "\t")
}


# capture session info
log("\n", func = "cli_text")
log("Session information:", func = "cli_h2")
si <- capture.output(sessionInfo())
for (i in si) log(i, func = "cli_text")

