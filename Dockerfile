# linux base with R version 4.3.3
FROM --platform=linux/amd64 rocker/r-ver:4.3.3

# linux installs
RUN apt-get update && \
    apt-get install -y \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    build-essential \
    libclang-dev \
    libcurl4-openssl-dev \
    libomp-dev \
    python3 \
    python3-pip \
    graphviz \
    cmake \
    libpng-dev \
    wget \
    unzip \
    make \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*


## pip installs
RUN pip3 install snakemake==7.26
RUN pip3 install --force-reinstall pulp==2.7.0

# plink2 install
WORKDIR /opt
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20260110.zip && \
    mkdir PLINK2 && \
    unzip plink2_linux_amd_avx2_20260110.zip -d PLINK2 && \
    rm plink2_linux_amd_avx2_20260110.zip
ENV PATH="/opt/PLINK2:${PATH}"

# R installs
RUN R -e "install.packages('data.table', version='1.17.8', type = 'source')"

RUN R -e 'install.packages(c( \
  "R.utils", \
  "viridis", \
  "ggplot2", \
  "ggpubr", \
  "GGally", \
  "rmarkdown", \
  "knitr", \
  "kableExtra", \
  "plotly", \
  "htmlwidgets" \
), repos="https://cloud.r-project.org")'

RUN R -e 'install.packages(c( \
  "lmerTest", \
  "partR2", \
  "rsq", \
  "confintr", \
  "future", \
  "furrr", \
  "progressr", \
  "broom", \
  "broom.mixed" \
), repos="https://cloud.r-project.org")'

RUN R -e 'install.packages(c( \
  "pheatmap", \
  "dynamicTreeCut", \
  "NbClust", \
  "cluster", \
  "fpc", \
  "circlize", \
  "eulerr", \
  "ggvenn" \
), repos="https://cloud.r-project.org")'

RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); \
           BiocManager::install(c( \
             "Rsamtools", \
             "SummarizedExperiment", \
             "MatrixGenerics", \
             "DelayedArray", \
             "GenomicAlignments", \
             "rtracklayer" \
           ), ask=FALSE)'

RUN R -e 'options( \
    repos = c(universe = "https://mrcieu.r-universe.dev/bin/linux/jammy/4.3/", \
              CRAN     = "https://packagemanager.posit.co/cran/__linux__/jammy/latest") \
  ); \
  install.packages("TwoSampleMR", dependencies = TRUE); \
  install.packages("genepi.utils", dependencies = TRUE)'

# other R installs
RUN R -e 'install.packages("readxl", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("lubridate", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("cowplot", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("ggvenn", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("RColorBrewer", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("GenSA")'
RUN R -e 'install.packages("eulerr", repos="https://cloud.r-project.org")'




# run this when finished creating
CMD ["/bin/bash"]

# build
# cd this_Dockerfile/dir/path
# docker build --progress=plain --tag nicksunderland/bmi_metabolites_and_heart_failure --file Dockerfile . 2>&1 | tee -a docker_build.log

# push to dockerhub
# docker push nicksunderland/bmi_metabolites_and_heart_failure:latest
