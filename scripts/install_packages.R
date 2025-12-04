# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install remotes if not available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# List of CRAN packages
cran_packages <- c(
  "devtools",
  "data.table",
  "R.utils",
  "lmerTest",
  "partR2",
  "rsq",
  "future",
  "furrr",
  "progressr",
  "broom.mixed",
  "broom",
  "viridis",
  "kableExtra",
  "rmarkdown",
  "pandoc",
  "htmlwidgets",
  "ggplot2",
  "GGally",
  "ggpubr",
  "ggvenn",
  "plotly",
  "reticulate",
  "nFactors",
  "pheatmap",
  "dynamicTreeCut",
  "NbClust",
  "cluster",
  "fpc",
  "multiUS",
  "mcr",
  "confintr",
  "pbmcapply",
  "circlize",
  "eulerr"
)


# Install missing CRAN packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
invisible(lapply(cran_packages, install_if_missing))


# stuff for MetaboAnalystR
install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/qvalue_2.40.0.tar.gz", type = "source", repos = NULL)
install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/SSPA_2.22.1.tar.gz", type = "source", repos = NULL)
metr_pkgs <- c(
    "impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz",
    "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph",
    "siggenes", "BiocParallel", "MSnbase", "multtest", "RBGL",
    "edgeR", "fgsea", "crmn", "edgeR"
)
installed    <- rownames(installed.packages())
missing_pkgs <- setdiff(metr_pkgs, installed)

if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  message("Installation complete.")
} else {
  message("No new packages needed.")
}

# GitHub packages, incl. MetaboAnalystR
github_packages <- c(
  "nicksunderland/metaboprep2",
  "xia-lab/MetaboAnalystR"
)

# Install GitHub packages
invisible(lapply(github_packages, remotes::install_github))
