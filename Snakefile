"""
Project: BMI metabolomics in heart failure
Author: nicholas.sunderland@bristol.ac.uk
Date: March 2025
"""

"""
SECTION: Global configuration and environment setup
"""
singularity: "docker://nicksunderland/bmi_metabolites_and_heart_failure:latest"
import os, glob
from snakemake.io import expand
from dotenv import load_dotenv
load_dotenv() # load the .env file variables - access environment variables with os.environ.get("VAR_NAME")
configfile: "data_config.yaml"
# resolve the GWAS catalog IDs
def resolve_gwas_paths(gwas_cfg):
    for name, cfg in gwas_cfg.items():
        if name == "default":
            continue

        # get base dir from environment variable
        dir_env_var = cfg.get("dir_env_var")
        if dir_env_var and not os.environ.get(dir_env_var):
            raise ValueError(f"Environment variable '{dir_env_var}' not set for '{name}'")
        base_dir = os.environ.get(dir_env_var, "") if dir_env_var else ""

        # resolve paths and ids
        id_builder = cfg.get("id_builder")
        if id_builder == "None":
            id_builder = None

        if id_builder == "glob":
            if not os.path.isdir(base_dir):
                raise ValueError(f"Directory does not exist for '{name}': '{base_dir}'")
            found = sorted(glob.glob(os.path.join(base_dir, "*.fst")))
            if not found:
                raise ValueError(f"No .fst files found in '{base_dir}' for '{name}'")
            cfg["paths"] = found
            cfg["ids"]   = [os.path.basename(f).replace(".fst", "") for f in found]

        elif id_builder is None:
            raw = cfg.get("paths")
            raw = raw if isinstance(raw, list) else []
            if not raw:
                raise ValueError(f"No paths defined for '{name}'")
            cfg["paths"] = [os.path.join(base_dir, p) for p in raw]
            cfg["ids"]   = [name]
            for p in cfg["paths"]:
                if not os.path.exists(p):
                    raise ValueError(f"Path does not exist for '{name}': '{p}'")
resolve_gwas_paths(config["gwases"])

def resolve_all_gwases(gwas_cfg):
    """Dict of {id: {platform, id, path, build}} across all gwases."""
    all_gwases = {}
    for name, cfg in gwas_cfg.items():
        if name == "default":
            continue
        paths = cfg.get("paths") or []
        ids   = cfg.get("ids")   or []
        build = cfg.get("gwas_build", "unknown")
        for i, p in zip(ids, paths):
            all_gwases[i] = {"platform": name, "id": i, "path": p, "build": build}
    return all_gwases
ALL_GWASES  = resolve_all_gwases(config["gwases"])
Hg19_GWASES = {i: e for i, e in ALL_GWASES.items() if e["build"] == "Hg19"}
Hg38_GWASES = {i: e for i, e in ALL_GWASES.items() if e["build"] == "Hg38"}
Hg19_BASE   = next(iter(Hg19_GWASES.values())) if Hg19_GWASES else None
Hg38_BASE   = next(iter(Hg38_GWASES.values())) if Hg38_GWASES else None
BMI_GWAS          = {i: e for i, e in ALL_GWASES.items() if i == "bmi"}
OUTCOME_GWASES    = {i: e for i, e in ALL_GWASES.items() if i in ["heart_failure","hip_knee_osteoarthritis","endometrial_cancer"]}
METABOLITE_GWASES = {i: e for i, e in ALL_GWASES.items() if e["platform"] in ["ms_chen_metabs","nmr_karjalainen_metabs"]}

# enumerate valid (pathway, gwas_name, gene_name) triples from use_with
PATHWAY_CIS_COMBOS = [
    (pathway.strip(), gwas_name.strip(), gene_name.strip())
    for pathway, p_data in config["pathway_instrument_search"].items()
    for gwas_name, g_data in p_data["gwases"].items()
    for gene_name in g_data.get("use_with", [])
    if gene_name.strip() in p_data["genes"]
]
# pathways that have at least one valid (gwas, gene) combo
PATHWAYS_WITH_COMBOS = list(dict.fromkeys(p for p, g, gn in PATHWAY_CIS_COMBOS))


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
SECTION: Cross-study replication and comparison
"""
checkpoint study_replication:
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
    name_map     = "scripts/gwas_metab_name_map.xlsx",
  output:
    assoc_table    = os.path.join("output", "tables", "linear_mixed_associations", "combined_trial_associations.tsv"),
    vol_bbs        = os.path.join("output", "figures", "associations", "bbs_linear_mixed_volcanos.png"),
    vol_direct     = os.path.join("output", "figures", "associations", "direct_linear_mixed_volcanos.png"),
    plat_corr      = os.path.join("output", "figures", "replication",  "platform_measurement_correlation.png"),
    plat_bland_alt = os.path.join("output", "figures", "replication",  "platform_measurement_bland_alterman.png"),
    corr_plot      = os.path.join("output", "figures", "replication",  "study_comparison_correlation.png"),
    study_bland_alt= os.path.join("output", "figures", "replication",  "study_comparison_bland_alterman.png"),
    rep_plot       = os.path.join("output", "figures", "replication",  "study_comparison_replication.png"),
    rep_table           = os.path.join("output", "tables",  "replication",  "study_replication.tsv"),
    consistent_ids      = os.path.join("output", "tables",  "replication",  "consistent_gwas_ids.txt"),
    direct_est_scatter  = os.path.join("output", "figures", "replication",  "direct_estimate_comparison_scatter.png"),
    direct_est_bland_alt= os.path.join("output", "figures", "replication",  "direct_estimate_comparison_bland_altman.png"),
    direct_se_scatter   = os.path.join("output", "figures", "replication",  "direct_se_comparison_scatter.png")
  log:
    log = os.path.join("output", "logs", "study_replication.log")
  script:
    "scripts/01_study_replication.R"






"""
SECTION: VARIANT MAPPING
"""
rule make_variant_map_extract:
    resources:
        cpus_per_task = 1,
        mem_mb        = 2000,
        time          = "00:02:00"
    input:
        input_gwas_fp = lambda wc: ALL_GWASES[wc.id]["path"]
    params:
        input_gwas_build = lambda wc: ALL_GWASES[wc.id]["build"],
        input_gwas_map   = lambda wc: config["gwases"][ALL_GWASES[wc.id]["platform"]]["map"],
        base_hg19_fp     = Hg19_BASE["path"],
        base_hg38_fp     = Hg38_BASE["path"],
        base_hg19_map    = config["gwases"][Hg19_BASE["platform"]]["map"],
        base_hg38_map    = config["gwases"][Hg38_BASE["platform"]]["map"]
    output:
        out = "output/tmp_objects/mapping/{id}_variants.tsv.gz"
    script:
        "scripts/make_variant_map_extract_step.R"

rule make_variant_map_merge:
    resources:
        cpus_per_task = 6,
        mem_mb        = 100000,
        time          = "05:59:00"
    input:
        variant_files = expand("output/tmp_objects/mapping/{id}_variants.tsv.gz", id = list(ALL_GWASES.keys()))
    params:
        dbsnp_dir = os.environ.get("DBSNP_DIR")
    output:
        map_file = "output/tables/global_variant_map.tsv.gz"
    script:
        "scripts/make_variant_map_merge_step.R"

"""
SECTION: CREATE INSTRUMENTS
"""
rule make_genetic_instruments:
  resources:
    cpus_per_task = 2,
    mem_mb        = 5000,
    time          = "00:10:00"
  input:
    gwas_file    = lambda wc: ALL_GWASES[wc.id]["path"],
    mapping_file = "output/tables/global_variant_map.tsv.gz"
  params:
    trait    = lambda wc: wc.id,
    id       = lambda wc: os.path.basename(ALL_GWASES[wc.id]["path"]),
    map      = lambda wc: config["gwases"][ALL_GWASES[wc.id]["platform"]]["map"],
    platform = lambda wc: ALL_GWASES[wc.id]["platform"],
    build    = lambda wc: ALL_GWASES[wc.id]["build"],
    r2       = 0.001,
    p1       = 5e-8,
    kb       = 10000,
    pfile = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k"),
  output:
    clump_file = "output/tables/instruments/genome_wide/{id}_instrument.tsv"
  script:
    "scripts/clump_gwas.R"


"""
SECTION: INSTRUMENT CLUSTERING 
"""
def consistent_metabolite_instrument_inputs(wildcards):
    co = checkpoints.study_replication.get()
    ids = [i for i in open(co.output.consistent_ids).read().splitlines() if i]
    metab_ids = [i for i in ids if i in METABOLITE_GWASES]
    return expand("output/tables/instruments/genome_wide/{id}_instrument.tsv", id=metab_ids)

rule instrument_clusters:
  input:
    instrument_files    = consistent_metabolite_instrument_inputs,
    consistent_ids_file = os.path.join("output", "tables", "replication", "consistent_gwas_ids.txt"),
    name_map_file       = "scripts/gwas_metab_name_map.xlsx",
    ld_block_file       = "scripts/ld_blocks_hg19.tsv"
  output:
    overlap_fig     = os.path.join("output", "figures", "instruments", "overlap_graph.png"),
    heatmap_fig     = os.path.join("output", "figures", "instruments", "overlap_heatmap.png"),
    ari_fig         = os.path.join("output", "figures", "instruments", "cluster_ari_stability.png"),
    cluster_table   = os.path.join("output", "tables",  "instruments", "cluster_membership.tsv"),
    cluster_bar_fig = os.path.join("output", "figures", "instruments", "cluster_bar.png"),
    ld_blk_ari_fig      = os.path.join("output", "figures", "instruments", "ld_block_ari_stability.png"),
    ld_blk_overlap_fig  = os.path.join("output", "figures", "instruments", "ld_block_overlap_graph.png")
  script:
    "scripts/02_instrument_clusters.R"


"""
SECTION: RUN MENDELIAN RANDOMIZATION
"""
rule run_mr:
  resources:
    cpus_per_task = 2,
    mem_mb        = 15000,
    time          = "02:00:00"
  input:
    instrument_file = "output/tables/instruments/genome_wide/{exp_id}_instrument.tsv",
    outcome_file    = lambda wc: ALL_GWASES[wc.out_id]["path"],
    mapping_file    = "output/tables/global_variant_map.tsv.gz"
  params:
    exp_id         = lambda wc: wc.exp_id,
    exp_binary     = lambda wc: config["gwases"][ALL_GWASES[wc.exp_id]["platform"]]["binary"],
    exp_ncase      = lambda wc: config["gwases"][ALL_GWASES[wc.exp_id]["platform"]]["ncase"],
    exp_ncontrol   = lambda wc: config["gwases"][ALL_GWASES[wc.exp_id]["platform"]]["ncontrol"],
    exp_prevalence = lambda wc: config["gwases"][ALL_GWASES[wc.exp_id]["platform"]]["prevalence"],
    out_id         = lambda wc: wc.out_id,
    outcome_map    = lambda wc: config["gwases"][ALL_GWASES[wc.out_id]["platform"]]["map"],
    out_build      = lambda wc: ALL_GWASES[wc.out_id]["build"],
    out_binary     = lambda wc: config["gwases"][ALL_GWASES[wc.out_id]["platform"]]["binary"],
    out_ncase      = lambda wc: config["gwases"][ALL_GWASES[wc.out_id]["platform"]]["ncase"],
    out_ncontrol   = lambda wc: config["gwases"][ALL_GWASES[wc.out_id]["platform"]]["ncontrol"],
    out_prevalence = lambda wc: config["gwases"][ALL_GWASES[wc.out_id]["platform"]]["prevalence"],
  output:
    results_file = "output/tmp_objects/mr_results/{exp_id}_on_{out_id}_mr.tsv"
  script:
    "scripts/run_mr.R"

"""
SECTION: CONCATENATE MR RESULTS
"""
def consistent_bmi_metab_inputs(wildcards):
    co = checkpoints.study_replication.get()
    ids = [i for i in open(co.output.consistent_ids).read().splitlines() if i]
    return expand(
        "output/tmp_objects/mr_results/{exp}_on_{out}_mr.tsv",
        exp = list(BMI_GWAS.keys()),
        out = ids,
    )

rule concat_mr_bmi_to_metabolites:
  input: consistent_bmi_metab_inputs
  output: "output/tables/mr_results/bmi_on_metabolites.tsv"
  shell:
    r"""
    Rscript -e "library(data.table); f <- commandArgs(TRUE); rbindlist(lapply(f[-length(f)], fread), fill=TRUE, use.names=TRUE) |> fwrite(f[length(f)], sep='\t')" {input} {output}
    """

def consistent_metab_outcome_inputs(wildcards):
    co = checkpoints.study_replication.get()
    ids = [i for i in open(co.output.consistent_ids).read().splitlines() if i]
    return expand(
        "output/tmp_objects/mr_results/{exp}_on_{out}_mr.tsv",
        exp = ids,
        out = list(OUTCOME_GWASES.keys()),
    )

rule concat_mr_metabolites_to_outcomes:
  input: consistent_metab_outcome_inputs
  output: "output/tables/mr_results/metabolites_on_outcomes.tsv"
  shell:
    r"""
    Rscript -e "library(data.table); f <- commandArgs(TRUE); rbindlist(lapply(f[-length(f)], fread), fill=TRUE, use.names=TRUE) |> fwrite(f[length(f)], sep='\t')" {input} {output}
    """

def consistent_outcome_metab_inputs(wildcards):
    co = checkpoints.study_replication.get()
    ids = [i for i in open(co.output.consistent_ids).read().splitlines() if i]
    return expand(
        "output/tmp_objects/mr_results/{exp}_on_{out}_mr.tsv",
        exp = list(OUTCOME_GWASES.keys()),
        out = ids,
    )

rule concat_mr_outcomes_to_metabolites:
  input: consistent_outcome_metab_inputs
  output: "output/tables/mr_results/outcomes_on_metabolites.tsv"
  shell:
    r"""
    Rscript -e "library(data.table); f <- commandArgs(TRUE); rbindlist(lapply(f[-length(f)], fread), fill=TRUE, use.names=TRUE) |> fwrite(f[length(f)], sep='\t')" {input} {output}
    """


"""
SECTION: VISUALISATION
"""
rule plot_bmi_metab_mr:
  input:
    results_file     = "output/tables/mr_results/bmi_on_metabolites.tsv",
    name_map_file    = "scripts/gwas_metab_name_map.xlsx",
    replicating_file = os.path.join("output", "tables", "replication", "study_replication.tsv")
  output:
    circos  = os.path.join("output", "figures", "bmi_mr", "bmi_metab_mr_circos.png"),
    scatter = os.path.join("output", "figures", "bmi_mr", "bmi_metab_mr_scatter.png")
  log:
    log = os.path.join("output", "logs", "bmi_metab_mr_discordance.log")
  script:
    "scripts/plot_bmi_metab_mr.R"

rule plot_metab_outcome_mr:
  input:
    results_file     = "output/tables/mr_results/metabolites_on_outcomes.tsv",
    rev_mr_file      = "output/tables/mr_results/outcomes_on_metabolites.tsv",
    name_map_file    = "scripts/gwas_metab_name_map.xlsx",
    replicating_file = os.path.join("output", "tables", "replication", "study_replication.tsv"),
    cluster_file     = os.path.join("output", "tables", "instruments", "cluster_membership.tsv")
  output:
    heatmap = os.path.join("output", "figures", "outcomes", "metab_outcome_mr.png"),
    circos  = os.path.join("output", "figures", "outcomes", "metab_outcome_mr_circos.png")
  script:
    "scripts/plot_metab_outcome_mr.R"


"""
SECTION: Project-wide target for full run
"""
rule all:
  input:
    # SECTION: CONCATENATE MR RESULTS
    "output/tables/mr_results/bmi_on_metabolites.tsv",
    "output/tables/mr_results/metabolites_on_outcomes.tsv",
    "output/tables/mr_results/outcomes_on_metabolites.tsv",
    # SECTION: INSTRUMENT CLUSTERING
    os.path.join("output", "figures", "instruments", "overlap_graph.png"),
    os.path.join("output", "figures", "instruments", "overlap_heatmap.png"),
    os.path.join("output", "figures", "instruments", "cluster_ari_stability.png"),
    # SECTION: VISUALISATION
    os.path.join("output", "figures", "bmi_mr",  "bmi_metab_mr_circos.png"),
    os.path.join("output", "figures", "bmi_mr",  "bmi_metab_mr_scatter.png"),
    os.path.join("output", "figures", "outcomes", "metab_outcome_mr.png"),
    os.path.join("output", "figures", "outcomes", "metab_outcome_mr_circos.png"),

    #"output/tables/instruments/genome_wide/bmi_instrument.tsv"
    #"output/tmp_objects/mapping/GCST90501056_variants.tsv.gz"
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
    # "output/figures/outcomes/metabolite_outcome_mr.png",
    # "output/figures/outcomes/metabolite_outcome_mr_forest.png",
    # "output/figures/outcomes/asparagine_mr.png"
    # expand("output/figures/cis_mr/{pathway}_gene_loci.png", pathway=list(config["pathway_instrument_search"].keys())),
    # expand("output/figures/cis_mr/{pathway}_mr_forest.png", pathway=PATHWAYS_WITH_COMBOS),
    # "output/tables/mr_results/pathway_cis_mr_results.tsv.gz",
    # expand("output/figures/cis_mr/{pathway}_mr_forest.png",pathway=PATHWAYS_WITH_COMBOS),
    # "output/figures/cis_mr/cis_vs_gw_mr_scatter.png",
    # "output/figures/cis_mr/observational_comparison.png"


#
# """
# SECTION: Workflow DAG generation
# """
# rule workflow_dag:
#   output: "workflow_dag.png"
#   shell: "snakemake --rulegraph --verbose | dot -Tpng > workflow_dag.png"
#
#


#
# """
# SECTION: Pathway gene locus plots (metabolite GWAS × cis-gene regions)
# """
# PATHWAY_TSS_WINDOW_KB=556 # PMID: 35527238
# rule plot_pathway_gene_loci:
#   """
#   For each pathway in pathway_instrument_search, plot all metabolite GWASes
#   against all cis-gene regions defined for that pathway.
#   Trigger with: snakemake output/figures/cis_mr/{pathway}_gene_loci.png
#   or expand over all pathways in rule all.
#   """
#   resources:
#     cpus_per_task = 4,
#     mem_mb        = 20000,
#     time          = "00:30:00"
#   input:
#     metab_dir = os.environ.get("METAB_DIR")
#   params:
#     gwases = lambda w: config["pathway_instrument_search"][w.pathway]["gwases"],
#     genes  = lambda w: config["pathway_instrument_search"][w.pathway]["genes"],
#     window_kb = PATHWAY_TSS_WINDOW_KB
#   output:
#     locus_plot = "output/figures/cis_mr/{pathway}_gene_loci.png"
#   script:
#     "scripts/plot_metab_gene_loci.R"
#
#
#
# """
# SECTION: CREATE CIS-PATHWAY INSTRUMENTS FROM pathway_instrument_search
# """
# rule make_pathway_cis_instrument:
#   """
#   For each valid (pathway, gwas_name, gene_name) triple in PATHWAY_CIS_COMBOS,
#   subset the metabolite GWAS to the cis-gene window and clump to produce an instrument.
#   """
#   resources:
#     cpus_per_task = 8,
#     mem_mb        = 40000,
#     time          = "01:00:00"
#   input:
#     gwas_file = lambda w: os.path.join(
#       os.environ.get("METAB_DIR"),
#       config["pathway_instrument_search"][w.pathway]["gwases"][w.gwas_name]["path"]
#     )
#   params:
#     gwas_build = lambda w: config["pathway_instrument_search"][w.pathway]["gwases"][w.gwas_name]["build"],
#     gene_chr   = lambda w: config["pathway_instrument_search"][w.pathway]["genes"][w.gene_name]["chr"],
#     gene_start = lambda w: config["pathway_instrument_search"][w.pathway]["genes"][w.gene_name]["start"],
#     gene_end   = lambda w: config["pathway_instrument_search"][w.pathway]["genes"][w.gene_name]["end"],
#     gene_win_kb= PATHWAY_TSS_WINDOW_KB,
#     r2         = 0.001,
#     p1         = 5e-8,
#     kb         = 250,
#     map        = config["gwases"]["metabolites"]["map"],
#     fill_rsid  = "b37_dbsnp156",
#     cores      = 2,
#     dbsnp_dir  = os.environ.get("DBSNP_DIR"),
#     pfile      = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k"),
#     trait      = lambda w: (w.gwas_name + "__" + w.gene_name).strip(),
#     id         = lambda w: (w.gwas_name + "__" + w.gene_name).strip()
#   output:
#     out_file = "output/tables/instruments/pathway_cis/{pathway}/{gwas_name}__{gene_name}_instrument.tsv"
#   script:
#     "scripts/make_drug_target_instrument.R"
#
#
# rule run_pathway_cis_mr:
#   resources:
#     cpus_per_task = 1,
#     mem_mb        = 15000,
#     time          = "00:45:00"
#   input:
#     out_list = [
#       os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["heart_failure"]["path"]),
#       os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["hfref"]["path"]),
#       os.path.join(os.environ.get("DATA_DIR"), config["gwases"]["hfpef"]["path"])
#     ],
#     out_map  = "output/tables/outcome_sorted_rsid_map.fst",
#     exp_gwas = "output/tables/instruments/pathway_cis/{pathway}/{gwas_name}__{gene_name}_instrument.tsv"
#   params:
#     out_maps = {
#       "heart_failure": config["gwases"]["heart_failure"]["map"],
#       "hfref":         config["gwases"]["hfref"]["map"],
#       "hfpef":         config["gwases"]["hfpef"]["map"]
#     },
#     pathway   = "{pathway}",
#     gwas_name = "{gwas_name}",
#     gene_name = "{gene_name}",
#     pfile     = os.path.join(os.environ.get("DATA_DIR"), "genome_reference", "ukb_reference_genome", "uk10k")
#   output:
#     mr_result = "output/tmp_objects/pathway_cis_mr/{pathway}/{gwas_name}__{gene_name}_mr_results.tsv.gz"
#   script:
#     "scripts/run_metabolite_pathway_outcome_mr.R"
#
#
# rule concat_pathway_cis_mr_results:
#     input:
#       ["output/tmp_objects/pathway_cis_mr/" + p + "/" + g + "__" + gn + "_mr_results.tsv.gz"
#        for p, g, gn in PATHWAY_CIS_COMBOS]
#     output:
#       "output/tables/mr_results/pathway_cis_mr_results.tsv.gz"
#     shell:
#       """
#       TMPDIR=$(mktemp -d)
#       for f in {input}; do zcat "$f" > "$TMPDIR/$(basename $f .gz).txt"; done
#       awk 'NR==1 || FNR>1' "$TMPDIR"/*.txt | gzip > {output}
#       rm -r "$TMPDIR"
#       """
#
#
# """
# SECTION: Observational vs MR comparison for cis-MR metabolites
# """
# rule plot_observational_comparison:
#   """
#   Scatter plot: Liu et al. observational HR (x) vs genome-wide MR OR (y).
#   Only Liu-significant metabolites (FDR p < 0.05) with a matching GW-MR
#   result are shown. Faceted by HF outcome (All HF, HFrEF, HFpEF).
#   """
#   resources:
#     cpus_per_task = 1,
#     mem_mb        = 4000,
#     time          = "00:15:00"
#   input:
#     gw_mr       = "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz",
#     name_map    = "scripts/gwas_metab_name_map.xlsx",
#     liu_fp      = "scripts/circhf_circhf-2023-010896_supp2.xlsx",
#     julkunen_fp = "scripts/julkunen_2023_summary_statistics.csv",
#     rep         = "output/tables/replication/metabolite_replication_bbs_direct_mr.tsv"
#   output:
#     comparison_plot = "output/figures/cis_mr/observational_comparison.png"
#   script:
#     "scripts/observational_comparison.R"
#
#
# """
# SECTION: Scatter plot — genome-wide vs cis-MR effect estimates
# """
# rule plot_cis_vs_gw_mr_scatter:
#   """
#   Single scatter plot across all pathways: GW MR estimate (x) vs cis-MR estimate (y).
#   IVW only; colour by metabolite, shape by gene (1, 2, 3 ...).
#   """
#   resources:
#     cpus_per_task = 1,
#     mem_mb        = 4000,
#     time          = "00:15:00"
#   input:
#     cis_mr = "output/tables/mr_results/pathway_cis_mr_results.tsv.gz",
#     gw_mr  = "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz"
#   params:
#     gwas_ids = {
#       gwas_name: os.path.splitext(os.path.basename(
#         config["pathway_instrument_search"][pathway]["gwases"][gwas_name]["path"]
#       ))[0]
#       for pathway  in config["pathway_instrument_search"]
#       for gwas_name in config["pathway_instrument_search"][pathway]["gwases"]
#     }
#   output:
#     scatter_plot = "output/figures/cis_mr/cis_vs_gw_mr_scatter.png"
#   script:
#     "scripts/plot_cis_vs_gw_mr.R"
#
#
# """
# SECTION: Per-pathway cis-MR forest plots
# """
# rule plot_pathway_cis_mr_forest:
#   """
#   Forest plot of cis-MR results for all metabolite x gene instruments within
#   a pathway, faceted by outcome (All HF, HFrEF, HFpEF).
#   Genome-wide instrument results are shown at the top for comparison.
#   One plot per pathway with at least one valid instrument (PATHWAYS_WITH_COMBOS).
#   """
#   resources:
#     cpus_per_task = 1,
#     mem_mb        = 4000,
#     time          = "00:15:00"
#   input:
#     mr_files = lambda w: [
#       "output/tmp_objects/pathway_cis_mr/" + p + "/" + g + "__" + gn + "_mr_results.tsv.gz"
#       for p, g, gn in PATHWAY_CIS_COMBOS if p == w.pathway
#     ],
#     gw_mr = "output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz"
#   params:
#     pathway  = "{pathway}",
#     gwas_ids = lambda w: {
#       gwas_name: os.path.splitext(os.path.basename(
#         config["pathway_instrument_search"][w.pathway]["gwases"][gwas_name]["path"]
#       ))[0]
#       for gwas_name in config["pathway_instrument_search"][w.pathway]["gwases"]
#     }
#   output:
#     forest_plot = "output/figures/cis_mr/{pathway}_mr_forest.png"
#   script:
#     "scripts/plot_metab_pathway_mr.R"
