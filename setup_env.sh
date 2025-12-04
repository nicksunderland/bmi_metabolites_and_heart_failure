#!/bin/bash
set -e

# vars
ENV_NAME="wt1_wp1_036_bmi_hf_metabolomics"
ENV_FILE="environment.yml"
MINICONDA_INSTALLER="Miniconda3-latest"

# detect OS
OS="$(uname)"
ARCH="$(uname -m)"
echo "Detected OS: $OS"
echo "Detected Architecture: $ARCH"

# installer
get_miniconda_installer() {
  if [[ "$OS" == "Darwin" ]]; then
    if [[ "$ARCH" == "arm64" ]]; then
      INSTALLER="Miniconda3-latest-MacOSX-arm64.sh"
    else
      INSTALLER="Miniconda3-latest-MacOSX-x86_64.sh"
    fi
  elif [[ "$OS" == "Linux" ]]; then
    if [[ "$ARCH" == "aarch64" ]]; then
      INSTALLER="Miniconda3-latest-Linux-aarch64.sh"
    else
      INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
    fi
  else
    echo "Unsupported OS: $OS"
    exit 1
  fi

  curl -LO "https://repo.anaconda.com/miniconda/$INSTALLER"
  bash "$INSTALLER" -b -p "$INSTALL_DIR"
  rm "$INSTALLER"

  export PATH="$INSTALL_DIR/bin:$PATH"
  echo "Miniconda installed at $INSTALL_DIR"
}

# check conda installed already?
if ! command -v conda &> /dev/null; then
  install_miniconda
else
  echo "Conda found: $(which conda)"
fi

# all set up?
eval "$(conda shell.bash hook)"

# create or update conda env
if conda env list | grep -q "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists. Updating environment..."
  conda env update -f "$ENV_FILE" --name "$ENV_NAME" || {
    echo "Failed to update environment."
    exit 1
  }
else
  echo "Creating new Conda environment from $ENV_FILE..."
  conda env create -f "$ENV_FILE" || {
    echo "Failed to create environment."
    exit 1
  }
fi

# activate
echo "Activating Conda environment..."
conda activate "$ENV_NAME"

# install R packages
echo "Running R package installation..."
Rscript scripts/install_packages.R