## Metabolomic signature of weight loss and association with heart failure
###### Project ID: wt1_wp1_036_bmi_hf_metabolomics
###### Contact: nicholas.sunderland@bristol.ac.uk

### Introduction
This repository contains a reproducible Snakemake workflow for setting up and executing the analyses behind the project “Metabolomic signature of weight loss and association with heart failure.”
The pipeline integrates metabolomic data from weight-loss interventions (DiRECT and By-Band-Sleeve trials) with Mendelian randomisation analyses to identify metabolic signatures of BMI change and evaluate their links to heart failure.

This README focuses on getting the computational environment ready, configuring the project, and running the workflow end-to-end so that results can be reproduced on the Bristol HPC. 

---

### TL;DR
Make a conda environment containing `snakemake==7.26`, `pulp==2.7.0`, and `python-dotenv`, clone the repository, set up the `.env` file (below), then run the analyses via the snakemake pipeline:
``` 
conda env create -n snakemake
conda activate snakemake
conda install python=3.11 snakemake==7.26 pulp==2.7.0 python-dotenv

git clone https://github.com/nicksunderland/bmi_metabolites_and_heart_failure $WORK/projects

module load apptainer
conda activate snakemake
cd $WORK/projects/bmi_metabolites_and_heart_failure
snakemake --profile . all
```

## Requirements

* `singularity 3.8.3`  
* `snakemake 7.26` (i.e. not 8.x)  
* `pulp 2.7.0` (i.e. <2.8)  
* Docker image: `docker://nicksunderland/bmi_metabolites_and_heart_failure:latest`  

## Docker image
All the required software and packages have been installed into a docker image. The image specification is contained in the `Dockerfile` incase this needs to be rebuilt or adjusted.

## Environment variables (.env file)
Once you have your cloned repo, rename the `.env_example` file to `.env` and edit the environment variables. Here you need to specify the raw data files, unfortunately these are not publicly available, and you would need to either obtain permissions to use them (see data availability statement) or run this at Bristol with appropriate account permissions.

| Variable                      | Description                                                                                                                                               |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| `SLURM_ACCOUNT`               | Your SLURM account ID for the HPC.                                                                                                                        |
| `DATA_DIR`                    | Root data directory on the HPC (GWASes and genomic reference files).                                                                                      |
| `METAB_DIR`                   | Parent directory for metabolomics datasets.                                                                                                               |
| `NIGHTINGALE_METAB_DIR`       | Nightingale (Karjalainen et al.) metabolomics subdirectory of  `METAB_DIR`.                                                                               |
| `METABOLON_METAB_DIR`         | Metabolon (Chen et al.) metabolomics subdirectory of `METAB_DIR`.                                                                                         |
| `DBSNP_DIR`                   | dbSNP reference data directory for use with [genepi.utils](https://nicksunderland.github.io/genepi.utils/).                                               |
| `BBS_CLINICAL`                | By-Band-Sleeve clinical file (`03_sample_clinical_data_all_cols.csv`).                                                                                    |
| `BBS_WITHDRAWAL`              | BBS withdrawals file (`WoC_20230608`).                                                                                                                    |
| `BBS_METABOLON_RDATA`         | BBS Metabolon RData file (`ReportData.Rdata`).                                                                                                            |
| `BBS_METABOLON_LOG`           | BBS Metabolon log file (`bbs_only_2025_04_07_logfile.txt`).                                                                                               |
| `BBS_METABOLON_MANIFEST`      | BBS Metabolon manifest file (`2023Q2_nmr_annotated_manifest_2023-08-14.csv`) — dummy variable required for workflow.                                      |
| `BBS_NIGHTINGALE_RDATA`       | BBS Nightingale RData file (`ReportData.Rdata`).                                                                                                          |
| `BBS_NIGHTINGALE_LOG`         | BBS Nightingale log file (`bbs_primary_nmr_no_outliers_2025_06_05_logfile.txt`).                                                                          |
| `BBS_NIGHTINGALE_MANIFEST`    | BBS Nightingale manifest file (`2023Q2_nmr_annotated_manifest_2023-08-14.csv`).                                                                           |
| `DIRECT_METABOLON`            | DiRECT Metabolon mixed-model results, run on the Glasgow TRE (`linear_mixed_associations_rnt_direct_scaled_metabolon.tsv`).                               |
| `DIRECT_INT_ONLY_METABOLON`   | DiRECT Metabolon intervention-only results, run on the Glasgow TRE (`linear_mixed_associations_rnt_direct_scaled_intervention_only_metabolon.tsv`).       |
| `DIRECT_NIGHTINGALE`          | DiRECT Nightingale mixed-model results, run on the Glasgow TRE (`linear_mixed_associations_raw_z_direct_scaled_nightingale.tsv`).                         |
| `DIRECT_INT_ONLY_NIGHTINGALE` | DiRECT Nightingale intervention-only results, run on the Glasgow TRE (`linear_mixed_associations_raw_z_direct_scaled_intervention_only_nightingale.tsv`). |

## Running the analyses
The project can be run on the Bristol HPC. This project uses a [Snakemake](https://snakemake.readthedocs.io/en/v7.32.0/getting_started/installation.html) pipeline to track analyses 
and ensure scripts and outputs are kept up-to-date. To reproduce the project end-to-end, first activate the conda environment and then run the Snakemake pipeline from this repo's root directory.  

First, load singularity and snakemake. Snakemake 7.26 is installed in a conda environment that needs to be activated.  

The project (slurm) configuration profile is contained within the `config.yaml` file. The snakemake `--profile` option takes a directory path to the where the `config.yaml` file can be found. In this case simply `.` as we have cd into the top level project directory.  

To run the project, submitting each job to the compute nodes, run the `all` rule like so:

```{bash}
module load apptainer
conda activate snakemake
cd $WORK/projects/bmi_metabolites_and_heart_failure
snakemake --profile . all -n
```

### Workflow rule graph
This executes the following work flows:  
![Workflow DAG](workflow_dag.png)


#### Steps
| Step | Purpose                                                                                                                         | Key script(s)                                                                                   | Primary outputs                                                                                                        |
|------|---------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| 1    | Preprocess By-Band-Sleeve metabolomics (Metabolon, Nightingale), apply QC, and save cleaned objects plus QC reports.            | `scripts/read_bbs_data_qc_summary.R`                                                          | `output/tables/qc/*.tsv`, `output/tmp_objects/bbs_*.RDS`                                                                |
| 2    | Build model-ready data frames (raw, z-scored, rank-normalised) for each platform.                                               | `scripts/make_model_dataframe.R`                                                              | `output/tmp_objects/model_df_*_bbs_*.tsv`, `output/tmp_objects/make_model_df_*_log.txt`                                |
| 3    | Run linear mixed association models for BBS data.                                                                               | `scripts/linear_mixed_associations_bbs.R`                                                     | `output/tables/linear_mixed_associations/linear_mixed_associations_*_bbs_*.tsv`, corresponding log files               |
| 4    | Import DiRECT linear mixed association results generated on the Glasgow TRE.                                                    | `scripts/read_direct_data.R`                                                                  | `output/tables/linear_mixed_associations/linear_mixed_associations_*_direct_*.tsv`                                     |
| 5    | Build a BMI instrument and metabolite chr:pos→rsID maps for Metabolon and Nightingale GWAS.                                     | `scripts/clump_gwas.R`, `scripts/make_rsid_map0.R`, `scripts/make_rsid_map.R`                  | `output/tables/instruments/bmi_instrument.tsv`, `output/tables/*_metabolite_rsid_map.fst`                              |
| 6    | Run BMI→metabolite Mendelian randomisation across all metabolites and collate results.                                         | `scripts/run_metabolite_mr.R`                                                                 | `output/tmp_objects/metabolite_bmi_mr/*.tsv`, `output/tables/mr_results/mr_results_bmi_metabolites.tsv`               |
| 7    | Compare BBS vs DiRECT associations and summarise replication (including BMI MR overlaps).                                      | `scripts/study_comparisons.R`                                                                 | Combined association tables, replication tables, and replication/association figures under `output/figures/replication` |
| 8    | Create rsID maps for heart failure outcomes and sorted metabolite GWAS (Metabolon, Nightingale).                               | `scripts/make_outcome_sorted_rsid_map.R`, `scripts/make_metabolite_sorted_rsid_map.R`          | `output/tables/outcome_sorted_rsid_map.fst`, `output/tables/{platform}_metabolite_sorted_rsid_map.fst`                 |
| 9    | Run metabolite→heart failure outcome Mendelian randomisation and generate summary tables/figures.                              | `scripts/run_metabolite_outcome_mr.R`, `scripts/metabolite_outcome_mr.R`                       | `output/tmp_objects/metabolite_outcome_mr/*.tsv.gz`, `output/tables/mr_results/metabolite_outcome_mr_results.tsv.gz`, outcome MR figures/tables             |
| 10   | Build pathway/gene target instruments and test their effects on heart failure outcomes.                                         | `scripts/make_drug_target_instrument.R`, `scripts/run_metabolite_pathway_outcome_mr.R`         | `output/tables/instruments/{pathway_target}_instrument.tsv`, `output/tables/mr_results/metabolite_pathway_outcome_mr_results.tsv.gz`                        |
| 11   | Plot asparagine-focused Mendelian randomisation summary.                                                                        | `scripts/asparagine_mr.R`                                                                     | `output/figures/outcomes/asparagine_mr.png`, `output/figures/outcomes/asparagine_mr_all_methods.png`                   |



---
