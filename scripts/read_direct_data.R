#### DESCRIPTION #######################################
# Purpose of script:
# Reads the DiRECT final results files exported from the
# Glasgow TRE.
#
# Author: Nick Sunderland
#
# Date Created: 2025-06-26
#
# Email: nicholas.sunderland@bristol.ac.uk

#### SET INPUT #########################################
file      <- snakemake@input[["file"]]
int_only  <- snakemake@input[["int_only"]]
map       <- snakemake@input[["map"]]
plat      <- snakemake@params[["platform"]]
data_type <- snakemake@params[["data_type"]]
main_study<- snakemake@params[["study"]]
assoc     <- snakemake@output[["assoc"]]
########################################################



# testing
if (FALSE) {
  plat <- "metabolon"
  #plat <- "nightingale"
  main_study <- "direct"
  if (plat == "nightingale") {
    file     <- "/Volumes/MRC-IEU-research/projects/wt1/wp1/036/working/data/direct_nick_results/linear_mixed_associations_raw_z_direct_scaled_nightingale.tsv"
    int_only <- "/Volumes/MRC-IEU-research/projects/wt1/wp1/036/working/data/direct_nick_results/linear_mixed_associations_raw_z_direct_scaled_intervention_only_nightingale.tsv"
  } else {
    file     <- "/Volumes/MRC-IEU-research/projects/wt1/wp1/036/working/data/direct_nick_results/linear_mixed_associations_rnt_direct_scaled_metabolon.tsv"
    int_only <- "/Volumes/MRC-IEU-research/projects/wt1/wp1/036/working/data/direct_nick_results/linear_mixed_associations_rnt_direct_scaled_intervention_only_metabolon.tsv"
  }
  data_type <- ifelse(plat=="metabolon", "rnt", "raw_z")
  assoc <- paste0("/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/.results/tables/linear_mixed_associations_", data_type, "_direct_", plat, ".tsv")
  map   <- "/Users/xx20081/git/wt1_wp1_036_bmi_hf_metabolomics/scripts/gwas_metab_name_map.xlsx"

}

# requirements
library(data.table)
library(readxl)


# read
interaction <- fread(file)
time_intervention_only <- fread(int_only)


# get interaction term
interaction <- interaction[term=="timepointend:treatIntervention"][, list(
  feature_id     = feature_id,
  term           = term,
  p.value        = p.value,
  p.value_fdr    = p.adjust(p.value, method = "fdr"),
  p.value_holm   = p.value_holm
)]


# get the timepoint term for the intervention only analysis
formatted <- time_intervention_only[term=="timepointend"][, list(
  feature_id     = feature_id,
  raw_id         = NA_character_,
  pathway        = NA_character_,
  sub_pathway    = NA_character_,
  derived_feature= NA,
  missingness    = NA_real_,
  outlier_count  = NA_integer_,
  independent    = NA,
  independent_k  = NA_integer_,
  effect         = effect,
  group          = group,
  term           = term,
  estimate       = estimate,
  std.error      = std.error,
  statistic      = statistic,
  df             = df,
  p.value        = p.value,
  rsq_tot        = rsq_tot,
  rsq_fixed      = rsq_fixed,
  rsq_random     = rsq_random,
  adj_rsq_tot    = adj_rsq_tot,
  adj_rsq_fixed  = adj_rsq_fixed,
  adj_rsq_random = adj_rsq_random,
  n              = n,
  main_study     = main_study,
  data_type      = data_type,
  platform       = platform,
  model_name     = model_name,
  formula        = formula,
  error	         = error,
  p.value_holm   = p.value_holm,
  p.value_fdr    = p.adjust(p.value, method = "fdr")
)]

formatted[interaction, `:=`(p.value_interaction      = i.p.value,
                            p.value_interaction_fdr  = i.p.value_fdr,
                            p.value_interaction_holm = i.p.value_holm), on="feature_id"]


# format for this project
if (plat == "nightingale") {

  # fix naming issues with nightingale ids in DiRECT
  # DiRECT -> map
  feature_map <- c(
    metab_vldld   = "metab_vldlsize",
    metab_ldld    = "metab_ldlsize",
    metab_hdld    = "metab_hdlsize",
    metab_serumc  = "metab_totalc",
    metab_estc    = "metab_totalce",
    metab_freec   = "metab_totalfc",
    metab_serumtg = "metab_totaltg",
    metab_totpg   = "metab_phosphoglyc",
    metab_pc      = "metab_phosphatidylc",
    metab_sm      = "metab_sphingomyelins",
    metab_totcho  = "metab_cholines",
    metab_totfa   = "metab_totalfa",
    metab_unsat   = "metab_unsaturation",
    metab_faw3    = "metab_omega3",
    metab_faw6    = "metab_omega6",
    metab_lafa    = "metab_lapct",
    metab_faw3fa  = "metab_omega3pct",
    metab_faw6fa  = "metab_omega6pct",
    metab_pufafa  = "metab_pufapct",
    metab_mufafa  = "metab_mufapct",
    metab_sfafa   = "metab_sfapct",
    metab_glc     = "metab_glucose",
    metab_lac     = "metab_lactate",
    metab_pyr     = "metab_pyruvate",
    metab_cit     = "metab_citrate",
    metab_glol    = "metab_glycerol",
    metab_ace     = "metab_acetate",
    metab_bohbut  = "metab_bohbutyrate",
    metab_crea    = "metab_creatinine",
    metab_alb     = "metab_albumin",
    metab_gp      = "metab_glyca"
  )
  formatted[, feature_id := ifelse(feature_id %in% names(feature_map), feature_map[feature_id], feature_id)]

}


# read the map and add standard name and pathway group
name_map <- read_xlsx(map, sheet=1) |> as.data.table()
formatted[name_map, `:=`(label = i.label, pathway_group = i.pathway_group), on=c("feature_id"="id")]
setcolorder(formatted, c("label","pathway_group"), after = "feature_id")


# save
fwrite(formatted, assoc, sep="\t")



# fwrite(formatted[is.na(label), "feature_id"], "~/Desktop/missing_direct.tsv", sep="\t")
# missing <- formatted[is.na(label), "feature_id"]
# missing[, compid := sub("metab_compid_", "", feature_id)]
#
# cache_dir <- tools::R_user_dir("metaboprep2", which = "cache")
# dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
# cache_file <- file.path(cache_dir, "master_compound_db.qs")
# db = data.table::as.data.table(qs::qread(cache_file))
#
# # compid ->
# new_map <- list(
#   '53174' = 'PE(18:2(9Z,12Z)/18:2(9Z,12Z))',
#   '57341' = 'PE(16:0/16:0)',
#   '57375' = 'PI(16:0/16:0)',
#   '35103' = 'LysoPC(20:2(11Z,14Z)/0:0)',
#   '36812' = 'LysoPA(18:1(9Z)/0:0)',
#   '52683' = 'PC(16:1(9Z)/18:2(9Z,12Z))',
#   '52684' = 'PC(16:0/18:3(9Z,12Z,15Z))',
#   '54812' = 'PC(16:0/18:3(6Z,9Z,12Z))',
#   '31928' = '2-Methylbutyrylglycine',
#   '52432' = 'MG(16:1(9Z)/0:0/0:0)',
#   '52929' = '3,4-methylene-heptanoylcarnitine',
#   '52938' = '3-Hydroxy stearic acid',
#   '32349' = 'Imidazoleacetic acid',
#   '43830' = '7-oxolithocholic acid',
#   '27414' = 'beta-Sitosterol',
#   '40007' = 'N-Carboxyethyl-g-aminobutyric acid',
#   '55002' = 'diacylclycerol (12:0/18:1, 14:0/16:1, 16:0/14:1) [1]*',
#   '32417' = 'Docosatrienoic acid',
#   '57467' = 'docosatrienoate (22:3n6)*',
#   '2134' = 'FAD',
#   '37092' = 'gamma-Glutamyl-2-aminobutyrate',
#   '61884' = 'glucuronide of C10H18O2 (4)*',
#   '48857' = 'Glycerophosphoglycerol',
#   '62065' = 'glyco-alpha-muricholate',
#   '57457' = 'glycosyl ceramide (d16:1/24:1, d18:1/22:1)*',
#   '57711' = 'HWESASXX*',
#   '32426' = 'Stercobilinogen',
#   '52749' = 'Indoxyl glucuronide',
#   '54949' = 'linoleoyl-docosahexaenoyl-glyerol (18:2/22:6) [1]*',
#   '54963' = 'DG(18:2(9Z,12Z)/18:3(6Z,9Z,12Z)/0:0)',
#   '37059' = 'Malonylcarnitine',
#   '35133' = '2-Methylguanosine',
#   '52632' = 'DG(16:1(9Z)/18:1(9Z)/0:0)',
#   '47695' = 'X - 12117',
#   '47659' = 'X - 12329',
#   '47922' = 'X - 12543',
#   '46318' = 'X - 12714',
#   '47938' = 'X - 12731',
#   '47984' = 'X - 14082',
#   '46656' = 'X - 14662',
#   '47275' = 'X - 16343',
#   '46490' = 'X - 46490',
#   '47997' = 'X - 17346',
#   '47009' = 'X - 17354',
#   '48019' = 'X - 17655',
#   '48022' = 'X - 17677',
#   '48024' = 'X - 48024',
#   '48055' = 'X - 17765',
#   '48049' = 'X - 18935',
#   '48050' = 'X - 18938',
#   '46316' = 'X - 21295',
#   '46359' = 'X - 21315',
#   '47772' = 'X - 22508',
#   '47773' = 'X - 22509',
#   '47776' = 'X - 22512',
#   '52659' = 'X - 22918',
#   '48570' = 'X - 23046',
#   '48771' = 'X - 23157',
#   '48773' = 'X - 23159',
#   '48774' = 'X - 23160',
#   '48892' = 'X - 23196',
#   '48993' = 'X - 23294',
#   '48995' = 'X - 23295',
#   '49461' = 'X - 23585',
#   '52297' = 'X - 24293',
#   '52544' = 'X - 24348',
#   '52848' = 'X - 24527',
#   '52863' = 'X - 24542',
#   '53065' = 'X - 24637',
#   '53097' = 'X - 24669',
#   '53266' = 'X - 24747',
#   '57756' = 'X - 24972',
#   '62488' = 'X - 25259',
#   '62609' = 'X - 25316',
#   '62689' = 'X - 25396',
#   '62697' = 'X - 25404'
# )
# db_subset <- rbind(
#   db[name %in% unname(new_map)],
#   data.table(name = unlist(unname(new_map)[!unname(new_map) %in% db$name])),
#   fill = TRUE)
# new_map_df <- data.table(feature_id = paste0("metab_compid_", names(new_map)),
#                          name = unlist(unname(new_map)))
# db_subset[new_map_df, feature_id := i.feature_id, on="name"]
#
#
# fwrite(db_subset[, .(
#   id = feature_id,
#   platform = "metabolon",
#   label = name,
#   pathway_group = NA_character_,
#   pathway = main_class,
#   subpathway = sub_class,
#   bbs_name = NA_character_,
#   direct_name = name,
#   ukbb_olink_name = NA_character_,
#   gwas_name = NA_character_,
#   gwas_id = NA_character_,
#   bbs_pathway = NA_character_,
#   direct_pathway = NA_character_,
#   bbs_subpathway = NA_character_,
#   direct_subpathway	 = NA_character_,
#   gwas_pathway = NA_character_,
#   gwas_subpathway = NA_character_
# )], "~/Desktop/map_updates_direct.tsv", sep="\t")
