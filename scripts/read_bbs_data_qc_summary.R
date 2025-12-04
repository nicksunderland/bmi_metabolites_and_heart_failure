#### DESCRIPTION #######################################
# Purpose of script:
# Reads various metabolite data formats, creates a Metaboprep2
# S7 class object holding the data in a standardised format
# and saves this object as a .RDS file.
#
# Author: Nick Sunderland
#
# Date Created: 2025-04-08
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
rdata         <- snakemake@input[["rdata"]]
clinical      <- snakemake@input[["clinical"]]
manifest      <- snakemake@input[["manifest"]]
map           <- snakemake@input[["map"]]
withdrawal    <- snakemake@input[["withdrawal"]]
olink         <- snakemake@input[["olink"]]
qc_summary    <- snakemake@output[["qc_summary"]]
qc_data_obj   <- snakemake@output[["qc_data_obj"]]
exc_samples   <- snakemake@output[["exc_samples"]]
exc_features  <- snakemake@output[["exc_features"]]
########################################################


# requirements
library(data.table)
library(readxl)


# load the metaboprep QC output
load(rdata)
clinic      <- fread(clinical)
withdrawals <- data.table(study_id = readLines(withdrawal)[!grepl("^#", readLines(withdrawal))])


# data.tables
raw_data$sample_data  <- as.data.table(raw_data$sample_data)
raw_data$feature_data <- as.data.table(raw_data$feature_data)
qc_data$sample_data   <- as.data.table(qc_data$sample_data)
qc_data$feature_data  <- as.data.table(qc_data$feature_data)


# get the clinical data in the samples data
if (platform=="Nightingale") {

  # rename
  raw_data$sample_data  <- raw_data$sample_data[,  .(sample_id       = Sample.id,
                                                     missingness     = sample_missingness,
                                                     outlier_count   = outlier_count)]
  qc_data$sample_data   <- qc_data$sample_data[,   .(sample_id       = Sample.id,
                                                     missingness     = sample_missingness,
                                                     outlier_count   = outlier_count)]
  raw_data$feature_data <- raw_data$feature_data[, .(feature_id      = paste0("metab_", tolower(gsub("_","", metabolite))),
                                                     raw_id          = metabolite,
                                                     name            = as.character(`Biomarker name`),
                                                     pathway         = as.character(Group),
                                                     sub_pathway     = as.character(Subgroup),
                                                     derived_feature = derived_features=="yes",
                                                     missingness     = feature_missingness,
                                                     outlier_count   = outlier_count,
                                                     independent     = independent_features_binary==1,
                                                     independent_k   = k)]
  qc_data$feature_data  <- qc_data$feature_data[,  .(feature_id      = paste0("metab_", tolower(gsub("_","", metabolite))),
                                                     raw_id          = metabolite,
                                                     name            = as.character(`Biomarker name`),
                                                     pathway         = as.character(Group),
                                                     sub_pathway     = as.character(Subgroup),
                                                     derived_feature = derived_features=="yes",
                                                     missingness     = feature_missingness,
                                                     outlier_count   = outlier_count,
                                                     independent     = independent_features_binary==1,
                                                     independent_k   = k)]


  # get the NMR manifest and add client sample id to samples
  stopifnot("The NMR manifest file must be provided if reading Nightingale NMR prerelease" = !is.null(manifest))
  manif <- fread(manifest)
  manif[, `:=`(barcode2   = sub("\\{\\}", "_", barcode),
               timepoint2 = data.table::fcase(grepl("36 months", timepoint), "end",
                                              grepl("oneyear", timepoint),   "end",
                                              grepl("(?i)baseline", timepoint),  "Baseline"))]
  manif[, combined_id := paste0(studyId, "_", timepoint2)]
  clinic[, combined_id := paste0(studyId, "_", timepoint)]
  clinic[manif, sample_id := i.barcode2, on="combined_id"]

  raw_data$sample_data[clinic, client_sample_id := i.Client.Sample.ID., on="sample_id"]
  qc_data$sample_data[clinic, client_sample_id := i.Client.Sample.ID., on="sample_id"]

} else if (platform=="Metabolon") {

  # rename
  raw_data$sample_data  <- raw_data$sample_data[,  .(sample_id       = PARENT_SAMPLE_NAME,
                                                     client_sample_id= CLIENT_SAMPLE_ID,
                                                     missingness     = sample_missingness,
                                                     outlier_count   = outlier_count)]
  qc_data$sample_data   <- qc_data$sample_data[,   .(sample_id       = PARENT_SAMPLE_NAME,
                                                     client_sample_id= CLIENT_SAMPLE_ID,
                                                     missingness     = sample_missingness,
                                                     outlier_count   = outlier_count)]
  raw_data$feature_data <- raw_data$feature_data[, .(feature_id      = paste0("metab_", feature_names),
                                                     raw_id          = feature_names,
                                                     name            = as.character(CHEMICAL_NAME),
                                                     pathway         = as.character(SUPER_PATHWAY),
                                                     sub_pathway     = as.character(SUB_PATHWAY),
                                                     derived_feature = FALSE,
                                                     missingness     = feature_missingness,
                                                     outlier_count   = outlier_count,
                                                     independent     = independent_features_binary==1,
                                                     independent_k   = k)]
  qc_data$feature_data  <- qc_data$feature_data[,  .(feature_id      = paste0("metab_", feature_names),
                                                     raw_id          = feature_names,
                                                     name            = as.character(CHEMICAL_NAME),
                                                     pathway         = as.character(SUPER_PATHWAY),
                                                     sub_pathway     = as.character(SUB_PATHWAY),
                                                     derived_feature = FALSE,
                                                     missingness     = feature_missingness,
                                                     outlier_count   = outlier_count,
                                                     independent     = independent_features_binary==1,
                                                     independent_k   = k)]
}


# get and clean clinical data
clinical_cols <- c(
  study_id         = "studyId",
  timepoint        = "timepoint",
  client_sample_id = "Client.Sample.ID.",
  study            = "source_study",
  sex              = "sex",
  age              = "age",
  ethnicity        = "ethnicity",
  height           = "hgt_m",
  weight           = "wgt",
  weight_op        = "wgt_op_day",
  bmi              = "bmikgm2",
  op_date          = "opdate",
  visit_date       = "visitdate",
  baseline_to_op   = "baseline_to_op",
  op_to_end        = "op_to_end",
  baseline_to_end  = "baseline_to_end",
  site             = "site",
  time_in_freezer  = "timeinfreezer"
)
clinic <- fread(clinical, select = unname(clinical_cols), col.names = names(clinical_cols))


# print follow up
fu_time <- clinic[, .(study_id, timepoint, op_date, visit_date, bmi)]
fu_time <- dcast(fu_time, study_id ~ timepoint, value.var = c("op_date", "visit_date", "bmi"))
fu_time[, fu_months := lubridate::interval(visit_date_Baseline, visit_date_end) / lubridate::dmonths(1)]
summary(fu_time)


# join clinical data
clinic[, join_flag := TRUE]
raw_data$sample_data <- clinic[raw_data$sample_data, on="client_sample_id", nomatch = NA]
qc_data$sample_data  <- clinic[qc_data$sample_data,  on="client_sample_id", nomatch = NA]


# assess missingness of clinical data
if (any(is.na(raw_data$sample_data$join_flag))) warning(paste0(sum(is.na(raw_data$sample_data$join_flag)), " raw samples with missing clinical data"))
if (any(is.na(qc_data$sample_data$join_flag))) warning(paste0(sum(is.na(qc_data$sample_data$join_flag)), " QC samples with missing clinical data"))


# add clinical data flags
clinic[, excl_withdrawal := study_id %in% withdrawals$study_id]
clinic[, excl_no_data := !study_id %in% raw_data$sample_data$study_id]
clinic[, excl_no_age := is.na(age)]
clinic[, excl_no_sex := is.na(sex)]


# read the map and add standard name and pathway group
name_map <- read_xlsx(map, sheet=1) |> as.data.table()
raw_data$feature_data[name_map, `:=`(label = i.label, pathway_group = i.pathway_group), on=c("feature_id"="id")]
qc_data$feature_data [name_map, `:=`(label = i.label, pathway_group = i.pathway_group), on=c("feature_id"="id")]
setcolorder(raw_data$feature_data, c("label","pathway_group"), after = "feature_id")
setcolorder(qc_data$feature_data, c("label","pathway_group"), after = "feature_id")


##################
# Quality control
##################

# individual exclusions
individual_exc_summary <- clinic[, .(
  individuals_total_n = uniqueN(study_id),
  individuals_main_total_n = uniqueN(study_id[study=="bbs_mainstudy"]),
  individuals_main_total_n = uniqueN(study_id[study=="bbs_substudy"]),
  withdrawal      = uniqueN(withdrawals$study_id),
  no_data         = uniqueN(study_id[excl_no_data==TRUE]),
  no_age          = uniqueN(study_id[excl_no_age==TRUE]),
  no_sex          = uniqueN(study_id[excl_no_sex==TRUE]),
  included        = uniqueN(study_id[excl_withdrawal==FALSE & excl_no_data==FALSE & excl_no_age==FALSE & excl_no_sex==FALSE])
)]
t(individual_exc_summary)


# sample exclusions - individual excluded due to withdrawal
raw_data$sample_data[, excl_withdrawal := study_id %in% withdrawals$study_id]
qc_data$sample_data[, excl_withdrawal := study_id %in% withdrawals$study_id]


# sample exclusions - individual excluded due to no sex data
raw_data$sample_data[, excl_no_sex := is.na(sex)]
qc_data$sample_data[, excl_no_sex := is.na(sex)]


# sample exclusions - individual excluded due to no age data
raw_data$sample_data[, excl_no_age := is.na(age)]
qc_data$sample_data[, excl_no_age := is.na(age)]


# sample exclusions - no operation
raw_data$sample_data[, excl_end_no_op := timepoint=="end" & is.na(op_date)]
qc_data$sample_data[, excl_end_no_op := timepoint=="end" & is.na(op_date)]


# sample exclusions - operation before end
raw_data$sample_data[, excl_end_b4_op := ifelse(is.na(op_to_end), FALSE, ifelse(timepoint=="end" & op_to_end <= 0, TRUE, FALSE))]
qc_data$sample_data[, excl_end_b4_op := ifelse(is.na(op_to_end), FALSE, ifelse(timepoint=="end" & op_to_end <= 0, TRUE, FALSE))]


# sample exclusion summary
sample_exc_summary <- cbind(
  qc_data$sample_data[, .(total_samps    = uniqueN(sample_id[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE]),
                          end_no_op      = uniqueN(sample_id[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==TRUE]),
                          end_b4_op      = uniqueN(sample_id[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_b4_op==TRUE]),
                          valid_n        = uniqueN(study_id[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==F & excl_end_b4_op==F]),
                          valid_samples  = uniqueN(sample_id[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==F & excl_end_b4_op==F]),
                          valid_baseline = uniqueN(sample_id[timepoint=="Baseline" & excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==F & excl_end_b4_op==F]),
                          valid_end      = uniqueN(sample_id[timepoint=="end" & excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==F & excl_end_b4_op==F])
                      )],
  data.table(valid_paired = sum(qc_data$sample_data[excl_withdrawal==FALSE & excl_no_sex==FALSE & excl_no_age==FALSE & excl_end_no_op==F & excl_end_b4_op==F,
                                     .(any(timepoint=="Baseline") & any(timepoint=="end")), by="study_id"]$V1))
)
t(sample_exc_summary)


# excluded sample ids
print(qcing_data$exclusion_data)
excluded_samples <- qc_data$sample_data[excl_withdrawal==TRUE | excl_no_sex==TRUE | excl_no_age==TRUE | excl_end_no_op==TRUE | excl_end_b4_op==TRUE]
excluded_samples
fwrite(excluded_samples, exc_samples, sep="\t")


# features
feature_miss_thresh <- 0.2
feature_exc_summary <- qc_data$feature_data[, .(features_total_n = uniqueN(feature_id),
                                                missing_gt20pct_n = uniqueN(feature_id[missingness>feature_miss_thresh]),
                                                valid_features_n  = uniqueN(feature_id[missingness<=feature_miss_thresh]))]
t(feature_exc_summary)


# feature ids
excluded_features <- qc_data$feature_data[missingness > feature_miss_thresh, ]
excluded_features
fwrite(excluded_features, exc_features, sep="\t")


# write QC summary
summ <- as.data.table(
  as.data.frame(
    t(cbind(individual_exc_summary, sample_exc_summary, feature_exc_summary))
  ),
  keep.rownames = TRUE)
setnames(summ, c("item", "count"))
summ
fwrite(summ, qc_summary, sep="\t")


###################
# Data to work with
###################

# extract metaboprep data
data     <- qc_data$metabolite_data
samples  <- qc_data$sample_data |> as.data.table()
features <- qc_data$feature_data |> as.data.table()


# exclude
samples  <- samples[!sample_id %in% excluded_samples$sample_id]
features <- features[!feature_id %in% excluded_features$feature_id]
data     <- data[which(rownames(data) %in% samples$sample_id), which(colnames(data) %in% features$raw_id)]
samples  <- samples[match(sample_id, rownames(data))]
features <- features[match(raw_id, colnames(data))]


# checks and return
stopifnot("sample ids don't match data rownames"  = identical(samples$sample_id, rownames(data)))
stopifnot("feature ids don't match data colnames" = identical(features$raw_id, colnames(data)))


# rename raw feature ids
colnames(data) <- features$feature_id


# create return list
ret_list <- list(
  data     = data,
  samples  = samples,
  features = features
)


# save the data object
saveRDS(ret_list, qc_data_obj)
