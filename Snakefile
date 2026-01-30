"""
Project: BMI metabolomics in heart failure
Author: nicholas.sunderland@bristol.ac.uk
Date: March 2025
"""

"""
SECTION: Global configuration and environment setup
"""
singularity: "docker://nicksunderland/bmi_metabolites_and_heart_failure:latest"
import os, sys
from snakemake.io import expand
from dotenv import load_dotenv
load_dotenv() # load the .env file variables - access environment variables with os.environ.get("VAR_NAME")
configfile: "data_config.yaml"
# load the metabolite GWAS catalog IDs
config["metabolite_ids"]["metabolon"] = ["GCST" + str(i) for i in range(90199621,90201021)]
config["metabolite_ids"]["nightingale"] = ["GCST" + str(i) for i in range(90301941,90302174)]


"""
SECTION: Project-wide target for full run
"""
rule all:
  input:
    # "output/tmp_objects/bbs_metabolon.RDS",
    # "output/tmp_objects/bbs_nightingale.RDS",
    # "output/tables/qc/bbs_nightingale_qc_summary.tsv",
    # "output/tables/qc/bbs_metabolon_qc_summary.tsv",
    # "output/tables/linear_mixed_associations/linear_mixed_associations_raw_z_bbs_nightingale.tsv",
    # "output/tables/linear_mixed_associations/linear_mixed_associations_rnt_bbs_metabolon.tsv",
    # "output/tables/linear_mixed_associations/linear_mixed_associations_raw_z_direct_nightingale.tsv",
    # "output/tables/linear_mixed_associations/linear_mixed_associations_rnt_direct_metabolon.tsv",
    # "output/figures/replication/study_comparison_bland_alterman.png",
    #"output/tmp_objects/mapping/metabolon_GCST90200194_rsid_map.tsv.gz"
    "output/figures/outcomes/metabolite_outcome_mr.png",
    "output/figures/outcomes/metabolite_outcome_mr_forest.png",
    "output/figures/outcomes/asparagine_mr.png"


"""
SECTION: Workflow DAG generation
"""
rule workflow_dag:
  output: "workflow_dag.png"
  shell: "snakemake --rulegraph --verbose | dot -Tpng > workflow_dag.png"


"""
SECTION: BBS data ingestion and QC summaries
"""
rule bbs_read_and_qc_summary:
  input:
    rdata    = lambda w: os.environ.get(f"BBS_{w.platform}".replace(" ", "").upper() + "_RDATA"),
    clinical = os.environ.get("BBS_CLINICAL"),
    manifest = lambda w: os.environ.get(f"BBS_{w.platform}".replace(" ", "").upper() + "_MANIFEST"),
    map      = "scripts/gwas_metab_name_map.xlsx",
    withdrawal = os.environ.get("BBS_WITHDRAWAL")
  output:
    qc_summary   = os.path.join("output", "tables", "qc", "{study}_{platform}_qc_summary.tsv"),
    exc_samples  = os.path.join("output", "tables", "qc", "{study}_{platform}_excluded_samples.tsv"),
    exc_features = os.path.join("output", "tables", "qc", "{study}_{platform}_excluded_features.tsv"),
    qc_data_obj  = os.path.join("output", "tmp_objects",  "{study}_{platform}.RDS")
  script:
    "scripts/read_bbs_data_qc_summary.R"


"""
SECTION: DiRECT mixed-model input ingestion
"""
rule read_direct_lmm:
  input:
    file      = lambda w: os.environ.get(f"DIRECT_{w.platform}".replace(" ", "").upper()),
    int_only  = lambda w: os.environ.get(f"DIRECT_INT_ONLY_{w.platform}".replace(" ", "").upper()),
    map       = "scripts/gwas_metab_name_map.xlsx",
  params:
    platform  = "{platform}",
    data_type = "{data_type}",
    study     = "direct"
  output:
    assoc     = os.path.join("output", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_direct_{platform}.tsv")
  script:
    "scripts/read_direct_data.R"


"""
SECTION: Model dataframe construction
"""
rule make_model_dataframes:
  input:
    qc_data_obj   = os.path.join("output", "tmp_objects", "{study}_{platform}.RDS")
  params:
    functions = "scripts/functions.R",
    platform  = "{platform}",
    study     = "{study}"
  output:
    model_raw          = os.path.join("output", "tmp_objects", "model_df_raw_{study}_{platform}.tsv"),
    model_raw_z        = os.path.join("output", "tmp_objects", "model_df_raw_z_{study}_{platform}.tsv"),
    model_rnt          = os.path.join("output", "tmp_objects", "model_df_rnt_{study}_{platform}.tsv"),
    log_file           = os.path.join("output", "tmp_objects", "make_model_df_{study}_{platform}_log.txt")
  script:
    "scripts/make_model_dataframe.R"


"""
SECTION: Linear mixed association models
"""
rule linear_mixed_associations:
  input:
    model_df      = os.path.join("output", "tmp_objects", "model_df_{data_type}_bbs_{platform}.tsv"),
    metaoprep_obj =  os.path.join("output", "tmp_objects", "bbs_{platform}.RDS")
  params:
    data_type = "{data_type}",
    platform  = "{platform}",
    functions = "scripts/functions.R"
  output:
    assoc    = os.path.join("output", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_bbs_{platform}.tsv"),
    log_file = os.path.join("output", "tables", "linear_mixed_associations", "linear_mixed_associations_{data_type}_bbs_{platform}_log.txt")
  script:
    "scripts/linear_mixed_associations_bbs.R"


"""
SECTION: CREATE THE BMI INSTRUMENT
"""
rule make_bmi_instrument:
  resources:
    cpus_per_task = 8,
    mem_mb        = 40000,
    time          = "01:00:00"
  input:
    gwas_file = os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["bmi"]["path"]),
  params:
    r2    = 0.001,
    p1    = 5e-8,
    kb    = 10000,
    map   = config["gwases"]["bmi"]["map"],
    pfile = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k"),
    trait = "bmi",
    id    = os.path.basename(config["gwases"]["bmi"]["path"])
  output:
    clump_file = "output/tables/instruments/bmi_instrument.tsv"
  script:
    "scripts/clump_gwas.R"


"""
SECTION: CREATE A MAPPING OF CHR:POS TO RSID FOR METABOLITE GWASES
"""
rule get_metabolite_gwas_variants:
  resources:
    cpus_per_task = 2,
    mem_mb        = 10000,
    time          = "00:30:00"
  input:
    base_file  = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ", "").upper()),
                                        config["metabolite_ids"][w.platform][0] + ".fst"),
    metab_file = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ", "").upper()),
                                        "{metab_gwas_id}.fst")
  output:
    out = "output/tmp_objects/mapping/{platform}_{metab_gwas_id}_rsid_map.tsv.gz"
  script:
    "scripts/make_rsid_map0.R"

rule make_rsid_map:
  resources:
    cpus_per_task = 12,
    mem_mb        = 150000,
    time          = "03:00:00"
  input:
    metab_files = lambda wildcards: expand("output/tmp_objects/mapping/{{platform}}_{id}_rsid_map.tsv.gz",
                                           id = config["metabolite_ids"][wildcards.platform]),
  params:
    dbsnp_dir   = os.environ.get("DBSNP_DIR"),
    platform    = "{platform}"
  output:
    map_file = "output/tables/{platform}_metabolite_rsid_map.fst"
  script:
    "scripts/make_rsid_map.R"


"""
SECTION: RUN MENDELIAN RANDOMIZATION BMI ON METABOLITE
"""
rule run_metabolite_mr:
  resources:
    cpus_per_task = 2,
    mem_mb = 20000,
    time = "00:15:00"
  input:
    exp_path = "output/tables/instruments/bmi_instrument.tsv",
    met_path = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ", "").upper()), f"{w.metabolite_id}.fst".replace(" ","")),
    met_yaml = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ", "").upper()), f"{w.metabolite_id}.tsv-meta.yaml".replace(" ","")),
    rsid_map = "output/tables/{platform}_metabolite_rsid_map.fst"
  params:
    platform       = "{platform}",
    metabolite_id  = "{metabolite_id}",
    metabolite_map = config["gwases"]["metabolites"]["map"],
    exposure       = "bmi",
    exp_clump_map  = config["gwases"]["default"]["map"],
    pfile          = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k")
  output:
    results_file = "output/tmp_objects/metabolite_bmi_mr/mr_results_bmi_{platform}_{metabolite_id}.tsv"
  script:
    "scripts/run_metabolite_mr.R"

rule concat_mr_results:
    input:
        [expand("output/tmp_objects/metabolite_bmi_mr/mr_results_bmi_{platform}_{id}.tsv",
                platform = "nightingale",
                id       = config["metabolite_ids"]["nightingale"]),
         expand("output/tmp_objects/metabolite_bmi_mr/mr_results_bmi_{platform}_{id}.tsv",
                platform = "metabolon",
                id       = config["metabolite_ids"]["metabolon"]),
         ]
    output:
        "output/tables/mr_results/mr_results_bmi_metabolites.tsv"
    shell:
        """
        awk 'FNR==1 && NR!=1 {{ next }} 1' {input} > {output}
        """


"""
SECTION: Cross-study replication and comparison
"""
rule study_comparisons:
  input:
    assoc = [expand("output/tables/linear_mixed_associations/linear_mixed_associations_{data_type}_{study}_{platform}.tsv",
                    data_type = ["raw_z"],
                    study     = ["bbs", "direct"],
                    platform  = ["nightingale"]),
             expand("output/tables/linear_mixed_associations/linear_mixed_associations_{data_type}_{study}_{platform}.tsv",
                    data_type = ["rnt"],
                    study     = ["bbs", "direct"],
                    platform = ["metabolon"])
             ],
    bbs_clinical = os.environ.get("BBS_CLINICAL"),
    name_map     = "scripts/gwas_metab_name_map.xlsx",
    bmi_mr       = "output/tables/mr_results/mr_results_bmi_metabolites.tsv"
  output:
    direct_dat     = os.path.join("output", "tables", "linear_mixed_associations", "combined_direct_associations.tsv"),
    bbs_dat        = os.path.join("output", "tables", "linear_mixed_associations", "combined_bbs_associations.tsv"),
    vol_bbs        = os.path.join("output", "figures", "associations", "bbs_linear_mixed_volcanos.png"),
    vol_direct     = os.path.join("output", "figures", "associations", "direct_linear_mixed_volcanos.png"),
    plat_corr      = os.path.join("output", "figures", "replication",  "platform_measurement_correlation.png"),
    plat_bland_alt = os.path.join("output", "figures", "replication",  "platform_measurement_bland_alterman.png"),
    corr_plot      = os.path.join("output", "figures", "replication",  "study_comparison_correlation.png"),
    study_bland_alt= os.path.join("output", "figures", "replication",  "study_comparison_bland_alterman.png"),
    rep_plot       = os.path.join("output", "figures", "replication",  "study_comparison_replication.png"),
    rep_table      = os.path.join("output", "tables",  "replication",  "metabolite_replication_bbs_direct_mr.tsv"),
    mr_rep_plot    = os.path.join("output", "figures", "replication",  "study_mr_replication.png"),
    circos_plot    = os.path.join("output", "figures", "replication",  "circos_plot.png")
  script:
    "scripts/study_comparisons.R"


"""
SECTION: RUN THE MENDELIAN RANDOMIZATION OF METABOLITES ON HEART FAILURE
"""
rule make_sorted_outcome_rsid_map:
  resources:
    cpus_per_task = 12,
    mem_mb        = 150000,
    time          = "01:00:00"
  input:
    gwases = [ os.path.join(os.environ["DATA_DIR"], config["gwases"]["heart_failure"]["path"]),
               os.path.join(os.environ["DATA_DIR"], config["gwases"]["hfref"]["path"]),
               os.path.join(os.environ["DATA_DIR"], config["gwases"]["hfpef"]["path"]) ]
  params:
    gwas_maps = {
      "heart_failure": config["gwases"]["heart_failure"]["map"],
      "hfref":         config["gwases"]["hfref"]["map"],
      "hfpef":         config["gwases"]["hfpef"]["map"],
    },
    dbsnp_dir = os.environ.get("DBSNP_DIR")
  output:
    map_out = "output/tables/outcome_sorted_rsid_map.fst"
  script:
    "scripts/make_outcome_sorted_rsid_map.R"

rule make_sorted_metabolite_rsid_map:
  resources:
    cpus_per_task = 12,
    mem_mb        = 150000,
    time          = "01:00:00"
  input:
    orig_map = "output/tables/{platform}_metabolite_rsid_map.fst"
  params:
    dbsnp_dir = os.environ.get("DBSNP_DIR")
  output:
    map_out = "output/tables/{platform}_metabolite_sorted_rsid_map.fst"
  script:
    "scripts/make_metabolite_sorted_rsid_map.R"

rule metabolite_outcome_mr:
  resources:
    cpus_per_task = 2,
    mem_mb        = 40000,
    time          = "03:00:00"
  input:
    out_list = [ os.path.join(os.environ["DATA_DIR"], config["gwases"]["heart_failure"]["path"]),
                 os.path.join(os.environ["DATA_DIR"], config["gwases"]["hfref"]["path"]),
                 os.path.join(os.environ["DATA_DIR"], config["gwases"]["hfpef"]["path"]) ],
    out_map     = "output/tables/outcome_sorted_rsid_map.fst",
    metab_gwas  = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ","").upper()),f"{w.metabolite_id}.fst".replace(" ","")),
    metab_meta  = lambda w: os.path.join(os.environ.get(f"{w.platform}_METAB_DIR".replace(" ","").upper()),f"{w.metabolite_id}.tsv-meta.yaml".replace(" ","")),
    metab_map   = "output/tables/{platform}_metabolite_sorted_rsid_map.fst"
  params:
    out_maps = {
      "heart_failure": config["gwases"]["heart_failure"]["map"],
      "hfref":         config["gwases"]["hfref"]["map"],
      "hfpef":         config["gwases"]["hfpef"]["map"],
    },
    platform = "{platform}",
    clump_r2 = 0.001,
    clump_p1 = 5e-8,
    clump_kb = 250,
    pfile    = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k")
  output:
    mr_instrument = "output/tmp_objects/metabolite_outcome_mr/{metabolite_id}_{platform}_outcome_mr_instrument.tsv.gz",
    mr_result     = "output/tmp_objects/metabolite_outcome_mr/{metabolite_id}_{platform}_outcome_mr_results.tsv.gz"
  script:
    "scripts/run_metabolite_outcome_mr.R"

rule concat_outcome_mr_results:
    input:
        [expand("output/tmp_objects/metabolite_outcome_mr/{id}_{platform}_outcome_mr_results.tsv.gz",
                platform = "nightingale",
                id       = config["metabolite_ids"]["nightingale"]),
         expand("output/tmp_objects/metabolite_outcome_mr/{id}_{platform}_outcome_mr_results.tsv.gz",
                platform = "metabolon",
                id       = config["metabolite_ids"]["metabolon"]),
         ]
    output:
        "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz"
    shell:
      """
      TMPDIR=$(mktemp -d)
      for f in {input}; do zcat "$f" > "$TMPDIR/$(basename $f .gz).txt"; done
      awk 'NR==1 || FNR>1' "$TMPDIR"/*.txt | gzip > {output}
      rm -r "$TMPDIR"
      """

rule concat_outcome_instruments:
    input:
        [expand("output/tmp_objects/metabolite_outcome_mr/{id}_{platform}_outcome_mr_instrument.tsv.gz",
                platform = "nightingale",
                id       = config["metabolite_ids"]["nightingale"]),
         expand("output/tmp_objects/metabolite_outcome_mr/{id}_{platform}_outcome_mr_instrument.tsv.gz",
                platform = "metabolon",
                id       = config["metabolite_ids"]["metabolon"]),
         ]
    output:
        "output/tables/mr_results/metabolite_outcome_instruments.tsv.gz"
    shell:
      """
      TMPDIR=$(mktemp -d)
      for f in {input}; do zcat "$f" > "$TMPDIR/$(basename $f .gz).txt"; done
      awk 'NR==1 || FNR>1' "$TMPDIR"/*.txt | gzip > {output}
      rm -r "$TMPDIR"
      """



"""
SECTION: Metabolite-outcome Mendelian randomisation
"""
rule plot_metabolite_outcome_mr:
  input:
    name_map      = "scripts/gwas_metab_name_map.xlsx",
    outcome_mr    = "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz",
    metab_loci_fp = "output/tables/mr_results/metabolite_outcome_instruments.tsv.gz",
    replicating   = os.path.join("output","tables","replication","metabolite_replication_bbs_direct_mr.tsv"),
    liu_fp        = "scripts/circhf_circhf-2023-010896_supp2.xlsx",
    julkunen_fp   = "scripts/julkunen_2023_summary_statistics.csv"
  params:
    outcome_snp_overlap_dir = os.path.join("output","figures","outcomes","metabolite_outcome_mr_sig_overlap")
  output:
    outcome_mr_tbl              = os.path.join("output", "tables", "outcomes", "metabolite_outcome_mr.tsv"),
    outcome_mr_plot             = os.path.join("output", "figures", "outcomes", "metabolite_outcome_mr.png"),
    outcome_mr_forest           = os.path.join("output", "figures", "outcomes", "metabolite_outcome_mr_forest.png"),
    outcome_volcano_forest      = os.path.join("output", "figures", "outcomes", "metabolite_outcome_volcano_forest.png"),
    outcome_mr_sig              = os.path.join("output", "figures", "outcomes", "metabolite_outcome_mr_sig.png"),
    outcome_instr_tbl           = os.path.join("output", "tables", "outcomes", "metabolite_outcome_instrument_overlap.tsv"),
    pathway_overlap_tbl         = os.path.join("output", "tables", "outcomes", "metabolite_outcome_instrument_pathway_overlap.tsv"),
    bmi_metabolite_snps         = os.path.join("output", "tables", "outcomes",  "metabolite_bmi_variants.tsv"),
    bmi_metabolite_variant_venn = os.path.join("output", "figures", "outcomes", "bmi_metabolite_variant_venn.png"),
    liu_model1_overlap_venn_ms  = os.path.join("output", "figures", "outcomes", "liu_model1_observational_metabolite_venn_ms.png"),
    liu_model1_overlap_venn_nmr = os.path.join("output","figures","outcomes","liu_model1_observational_metabolite_venn_nmr.png"),
    out_mr_obs_heatmap          = os.path.join("output", "figures", "outcomes", "outcomes/mr_observational_heatmap.png"),
    observational_hf_overlap_tbl= os.path.join("output", "tables", "outcomes", "observational_hf_overlap.tsv")
  script:
    "scripts/metabolite_outcome_mr.R"



"""
SECTION: CREATE A CIS-PATHWAY TARGET INSTRUMENT 
"""
rule make_gene_target_instrument:
  """
  This rule creates a pathway target instrument given parameters in the config file.
  """
  resources:
    cpus_per_task = 8,
    mem_mb        = 40000,
    time          = "01:00:00"
  input:
    gwas_file = lambda w: os.path.join(os.environ.get("METAB_DIR"), config["pathway_targets"][w.pathway_target]["gwas_path"])
  params:
    gwas_build = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["gwas_build"],
    gene_chr   = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["gene_chr"],
    gene_start = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["gene_start"],
    gene_end   = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["gene_end"],
    gene_win   = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["gene_win"],
    r2         = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["r2"],
    p1         = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["p1"],
    kb         = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["kb"],
    map        = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["map"],
    fill_rsid  = lambda wildcards: config["pathway_targets"][wildcards.pathway_target]["fill_rsid"],
    cores      = 2,
    dbsnp_dir  = os.environ.get("DBSNP_DIR"),
    pfile      = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k"),
    trait      = "{pathway_target}",
    id         = "{pathway_target}"
  output:
    out_file = "output/tables/instruments/{pathway_target}_instrument.tsv"
  script:
    "scripts/make_drug_target_instrument.R"


rule metabolite_gene_outcome_mr:
  resources:
    cpus_per_task = 1,
    mem_mb        = 15000,
    time          = "00:30:00"
  input:
    out_list = [
      os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["heart_failure"]["path"]),
      os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["hfref"]["path"]),
      os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["hfpef"]["path"])
    ],
    out_map     = "output/tables/outcome_sorted_rsid_map.fst",
    exp_gwas    = "output/tables/instruments/{pathway}_instrument.tsv"
  params:
    out_maps = {
        "heart_failure": config["gwases"]["heart_failure"]["map"],
        "hfref": config["gwases"]["hfref"]["map"],
        "hfpef": config["gwases"]["hfpef"]["map"]
    },
    pfile = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k"),
  output:
    mr_result = "output/tmp_objects/metabolite_outcome_mr/{pathway}_pathway_outcome_mr_results.tsv.gz"
  script:
    "scripts/run_metabolite_pathway_outcome_mr.R"


rule concat_gene_mr_results:
    input: expand("output/tmp_objects/metabolite_outcome_mr/{pathway}_pathway_outcome_mr_results.tsv.gz", pathway = list(config["pathway_targets"].keys()))
    output:
        "output/tables/mr_results/metabolite_pathway_outcome_mr_results.tsv.gz"
    shell:
      """
      TMPDIR=$(mktemp -d)
      for f in {input}; do zcat "$f" > "$TMPDIR/$(basename $f .gz).txt"; done
      awk 'NR==1 || FNR>1' "$TMPDIR"/*.txt | gzip > {output}
      rm -r "$TMPDIR"
      """


"""
SECTION: Asparagine Mendelian randomisation
"""
rule plot_gene_outcome_mr:
  input:
    metab_outcome_mr_results = "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz",
    metab_outcome_pathway_mr = "output/tables/mr_results/metabolite_pathway_outcome_mr_results.tsv.gz"
  output:
    asp_hf_mr_plot          = os.path.join("output", "figures", "outcomes", "asparagine_mr.png"),
    asp_hf_mr_all_plot      = os.path.join("output", "figures", "outcomes", "asparagine_mr_all_methods.png")
  script:
    "scripts/asparagine_mr.R"