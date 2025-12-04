"""
Project: BMI metabolomics in heart failure
Author: nicholas.sunderland@bristol.ac.uk
Date: March 2025
"""

"""
SECTION: Global configuration and environment setup
"""
container: "nicksunderland/wt1_wp1_036_bmi_hf_metabolomics"
import os, sys
from snakemake.io import expand
from dotenv import load_dotenv
load_dotenv() # load the .env file variables - access environment variables with os.environ.get("VAR_NAME")


"""
SECTION: Project-wide target for full run
"""
rule all:
  input:
    ".results/tmp_objects/bbs_metabolon.RDS",
    ".results/tmp_objects/bbs_nightingale.RDS",
    ".results/tables/qc/bbs_nightingale_qc_summary.tsv",
    ".results/tables/qc/bbs_metabolon_qc_summary.tsv",
    ".results/tables/linear_mixed_associations/linear_mixed_associations_raw_z_bbs_nightingale.tsv",
    ".results/tables/linear_mixed_associations/linear_mixed_associations_rnt_bbs_metabolon.tsv",
    ".results/tables/linear_mixed_associations/linear_mixed_associations_raw_z_direct_nightingale.tsv",
    ".results/tables/linear_mixed_associations/linear_mixed_associations_rnt_direct_metabolon.tsv",
    ".results/figures/replication/study_comparison_bland_alterman.png",
    ".results/figures/outcomes/metabolite_outcome_mr.png",
    ".results/figures/outcomes/metabolite_outcome_mr_forest.png",
    ".results/figures/outcomes/asparagine_mr.png"


"""
SECTION: Workflow DAG generation
"""
rule workflow_dag:
  """
  Workflow DAG generation: emit rule graph visual.
  Outputs:
  - workflow_dag.png: Graph of rule dependencies.
  """
  output: "workflow_dag.png"
  shell: "snakemake --rulegraph --verbose | dot -Tpng > workflow_dag.png"


"""
SECTION: BBS data ingestion and QC summaries
"""
rule bbs_read_and_qc_summary:
  input:
    rdata      = lambda wildcards: os.environ.get(f"BBS_{wildcards.platform.upper()}_RDATA"),
    clinical   = os.environ.get("BBS_CLINICAL"),
    manifest   = lambda wildcards: os.environ.get(f"BBS_{wildcards.platform.upper()}_MANIFEST"),
    map        = "scripts/gwas_metab_name_map.xlsx",
    withdrawal = os.environ.get("BBS_WITHDRAWAL")
  output:
    qc_summary   = os.path.join(".results", "tables", "qc", "{study}_{platform}_qc_summary.tsv"),
    exc_samples  = os.path.join(".results", "tables", "qc", "{study}_{platform}_excluded_samples.tsv"),
    exc_features = os.path.join(".results", "tables", "qc", "{study}_{platform}_excluded_features.tsv"),
    qc_data_obj  = os.path.join(".results", "tmp_objects",  "{study}_{platform}.RDS")
  script:
    "scripts/read_bbs_data_qc_summary.R"


"""
SECTION: DiRECT mixed-model input ingestion
"""
rule read_direct_lmm:
  input:
    file      = lambda wildcards: os.environ.get(f"DIRECT_{wildcards.platform.upper()}"),
    int_only  = lambda wildcards: os.environ.get(f"DIRECT_INT_ONLY_{wildcards.platform.upper()}"),
    map       = "scripts/gwas_metab_name_map.xlsx",
  params:
    platform  = "{platform}",
    data_type = "{data_type}",
    study     = "direct"
  output:
    assoc     = os.path.join(".results", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_direct_{platform}.tsv")
  script:
    "scripts/read_direct_data.R"


"""
SECTION: Model dataframe construction
"""
rule make_model_dataframes:
  input:
    qc_data_obj   = os.path.join(".results", "tmp_objects", "{study}_{platform}.RDS")
  params:
    functions = "scripts/functions.R",
    platform  = "{platform}",
    study     = "{study}"
  output:
    model_raw          = os.path.join(".results", "tmp_objects", "model_df_raw_{study}_{platform}.tsv"),
    model_raw_z        = os.path.join(".results", "tmp_objects", "model_df_raw_z_{study}_{platform}.tsv"),
    model_rnt          = os.path.join(".results", "tmp_objects", "model_df_rnt_{study}_{platform}.tsv"),
    log_file           = os.path.join(".results", "tmp_objects", "make_model_df_{study}_{platform}_log.txt")
  script:
    "scripts/make_model_dataframe.R"


"""
SECTION: Linear mixed association models
"""
rule linear_mixed_associations:
  input:
    model_df      = os.path.join(".results", "tmp_objects", "model_df_{data_type}_bbs_{platform}.tsv"),
    metaoprep_obj =  os.path.join(".results", "tmp_objects", "bbs_{platform}.RDS")
  params:
    data_type = "{data_type}",
    platform  = "{platform}",
    functions = "scripts/functions.R"
  output:
    assoc    = os.path.join(".results", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_bbs_{platform}.tsv"),
    log_file = os.path.join(".results", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_bbs_{platform}_log.txt")
  script:
    "scripts/linear_mixed_associations_bbs.R"



"""
SECTION: Cross-study replication and comparison
"""
rule study_comparisons:
  input:
    assoc = [expand(".results/tables/linear_mixed_associations/linear_mixed_associations_{data_type}_{study}_{platform}.tsv",
                    data_type = ["raw_z"],
                    study     = ["bbs", "direct"],
                    platform  = ["nightingale"]),
             expand(".results/tables/linear_mixed_associations/linear_mixed_associations_{data_type}_{study}_{platform}.tsv",
                    data_type = ["rnt"],
                    study     = ["bbs", "direct"],
                    platform = ["metabolon"])
             ],
    bbs_clinical = os.environ.get("BBS_CLINICAL"),
    name_map     = "scripts/gwas_metab_name_map.xlsx",
    bmi_mr       = os.path.join(os.environ.get("BMI_HF_MR_REPO"), "output", "tables", "mr_results", "mr_results_bmi_metabolites.tsv")
  output:
    direct_dat     = os.path.join(".results", "tables", "linear_mixed_associations", "combined_direct_associations.tsv"),
    bbs_dat        = os.path.join(".results", "tables", "linear_mixed_associations", "combined_bbs_associations.tsv"),
    vol_bbs        = os.path.join(".results", "figures", "associations", "bbs_linear_mixed_volcanos.png"),
    vol_direct     = os.path.join(".results", "figures", "associations", "direct_linear_mixed_volcanos.png"),
    plat_corr      = os.path.join(".results", "figures", "replication",  "platform_measurement_correlation.png"),
    plat_bland_alt = os.path.join(".results", "figures", "replication",  "platform_measurement_bland_alterman.png"),
    corr_plot      = os.path.join(".results", "figures", "replication",  "study_comparison_correlation.png"),
    study_bland_alt= os.path.join(".results", "figures", "replication",  "study_comparison_bland_alterman.png"),
    rep_plot       = os.path.join(".results", "figures", "replication",  "study_comparison_replication.png"),
    rep_table      = os.path.join(".results", "tables",  "replication",  "metabolite_replication_bbs_direct_mr.tsv"),
    mr_rep_plot    = os.path.join(".results", "figures", "replication",  "study_mr_replication.png"),
    circos_plot    = os.path.join(".results", "figures", "replication",  "circos_plot.png")
  script:
    "scripts/study_comparisons.R"


"""
SECTION: Metabolite-outcome Mendelian randomisation
"""
rule metabolite_outcome_mr:
  input:
    name_map      = "scripts/gwas_metab_name_map.xlsx",
    outcome_mr    = os.path.join(os.environ.get("BMI_HF_MR_REPO"), "output", "tables", "mr_results", "metabolite_outcome_mr_results.tsv.gz"),
    metab_loci_fp = os.path.join(os.environ.get("BMI_HF_MR_REPO"), "output", "tables", "mr_results", "metabolite_outcome_instruments.tsv.gz"),
    replicating   = os.path.join(".results","tables","replication","metabolite_replication_bbs_direct_mr.tsv")
  params:
    outcome_snp_overlap_dir = os.path.join(".results","figures","outcomes","metabolite_outcome_mr_sig_overlap")
  output:
    outcome_mr_tbl              = os.path.join(".results", "tables", "outcomes", "metabolite_outcome_mr.tsv"),
    outcome_mr_plot             = os.path.join(".results", "figures", "outcomes", "metabolite_outcome_mr.png"),
    outcome_mr_forest           = os.path.join(".results", "figures", "outcomes", "metabolite_outcome_mr_forest.png"),
    outcome_volcano_forest      = os.path.join(".results", "figures", "outcomes", "metabolite_outcome_volcano_forest.png"),
    outcome_mr_sig              = os.path.join(".results", "figures", "outcomes", "metabolite_outcome_mr_sig.png"),
    outcome_instr_tbl           = os.path.join(".results", "tables", "outcomes", "metabolite_outcome_instrument_overlap.tsv"),
    pathway_overlap_tbl         = os.path.join(".results", "tables", "outcomes", "metabolite_outcome_instrument_pathway_overlap.tsv"),
    bmi_metabolite_snps         = os.path.join(".results", "tables", "outcomes",  "metabolite_bmi_variants.tsv"),
    bmi_metabolite_variant_venn = os.path.join(".results", "figures", "outcomes", "bmi_metabolite_variant_venn.png"),
    liu_model1_overlap_venn_ms  = os.path.join(".results", "figures", "outcomes", "liu_model1_observational_metabolite_venn_ms.png"),
    liu_model1_overlap_venn_nmr = os.path.join(".results","figures","outcomes","liu_model1_observational_metabolite_venn_nmr.png"),
    out_mr_obs_heatmap          = os.path.join(".results", "figures", "outcomes", "outcomes/mr_observational_heatmap.png")
  script:
    "scripts/metabolite_outcome_mr.R"


"""
SECTION: Asparagine Mendelian randomisation
"""
rule asparagine_outcome_mr:
  input:
    metab_outcome_mr_results = os.path.join(os.environ.get("BMI_HF_MR_REPO"), "output", "tables", "mr_results", "metabolite_outcome_mr_results.tsv.gz"),
    metab_outcome_pathway_mr = os.path.join(os.environ.get("BMI_HF_MR_REPO"), "output", "tables", "mr_results", "metabolite_pathway_outcome_mr_results.tsv.gz")
  output:
    asp_hf_mr_plot          = os.path.join(".results", "figures", "outcomes", "asparagine_mr.png"),
    asp_hf_mr_all_plot      = os.path.join(".results", "figures", "outcomes", "asparagine_mr_all_methods.png")
  script:
    "scripts/asparagine_mr.R"