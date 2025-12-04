# linux base with R version 4.3.3
FROM --platform=linux/amd64 rocker/r-ver:4.3.3

# linux installs
RUN apt-get update && apt-get upgrade -y

RUN apt-get install -y \
    liblzma-dev \
    libbz2-dev \
    libclang-dev \
    libcurl4-openssl-dev \
    libomp-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    python3 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    python3-pip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    graphviz \
    cmake \
    gwama \
    wget \
    unzip \
    make \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*


## pip installs
RUN pip3 install snakemake==7.26
RUN pip3 install --force-reinstall pulp==2.7.0

# install VEP via conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

# Accept Anaconda Terms of Service for required channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create the environment and install VEP
RUN conda create -n vep -y && \
    conda install -n vep -c bioconda ensembl-vep -y

# Activate the environment by default
ENV PATH /opt/conda/envs/vep/bin:$PATH
ENV CONDA_DEFAULT_ENV vep

# plink2 install
WORKDIR /opt
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20251026.zip && \
    mkdir PLINK2 && \
    unzip plink2_linux_amd_avx2_20251026.zip -d PLINK2 && \
    rm plink2_linux_amd_avx2_20251026.zip
ENV PATH="/opt/PLINK2:${PATH}"

# R installs
RUN R -e "install.packages('fst', version='0.9.8', type = 'source')"
RUN R -e "install.packages('furrr', version='0.3.1', type = 'source')"
RUN R -e "install.packages('data.table', version='1.15.2', type = 'source')"
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('R.utils')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e 'options( \
    repos = c(universe = "https://mrcieu.r-universe.dev/bin/linux/jammy/4.3/", \
              CRAN     = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"), \
    HTTPUserAgent = sprintf( \
        "R/%s R (%s)", \
        getRversion(), \
        paste(getRversion(), \
          R.version["platform"], \
          R.version["arch"], \
          R.version["os"]))); \
    install.packages("TwoSampleMR", dependencies = TRUE); '
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('viridis')"
RUN R -e "install.packages('Boruta')"
RUN R -e "remotes::install_github('mglev1n/mrbma')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
           BiocManager::install(c('SummarizedExperiment', 'MatrixGenerics', 'DelayedArray', 'GenomicAlignments', 'rtracklayer'))"
RUN R -e "install.packages('caret')"
RUN R -e "install.packages('patchwork')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('rms')"
RUN R -e "install.packages('ggfortify')"
RUN R -e "install.packages('survminer')"
RUN R -e "install.packages('survivalAnalysis')"
RUN R -e "install.packages('finalfit')"
RUN R -e "install.packages('cmprsk')"
RUN R -e "install.packages('recipes')"
RUN R -e "remotes::install_github('nicksunderland/genepi.utils@v0.1.7')"

# LDSC
RUN apt-get update && apt-get install -y \
    git \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/bulik/ldsc.git /opt/ldsc
RUN conda env create -f /opt/ldsc/environment.yml
ENV PATH /opt/conda/envs/ldsc/bin:$PATH
ENV CONDA_DEFAULT_ENV ldsc
RUN /opt/ldsc/ldsc.py -h
RUN /opt/ldsc/munge_sumstats.py -h

# more packages
RUN R -e "install.packages('mets')"
RUN R -e "install.packages('prodlim')"
RUN R -e "install.packages('ggsurvfit')"
RUN R -e "install.packages('boot')"


# run this when finished creating
CMD ["/bin/bash"]

# build
# cd this_Dockerfile/dir/path
# docker build --progress=plain --tag nicksunderland/hf_mvmr --file Dockerfile . 2>&1 | tee -a docker_build.log

# push to dockerhub
# docker push nicksunderland/hf_mvmr:latest
