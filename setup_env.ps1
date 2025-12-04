# error
$ErrorActionPreference = "Stop"

# vars
$ENV_NAME = "wt1_wp1_036_bmi_hf_metabolomics"
$ENV_FILE = "environment.yml"

# Detect Windows Architecture
$ARCH = if ($env:PROCESSOR_ARCHITECTURE -eq "AMD64") { "x86_64" } else { "arm64" }
Write-Host "Detected Architecture: $ARCH"

# Function to Install Miniconda
function Install-Miniconda {
    Write-Host "Miniconda not found. Installing Miniconda..."

    # Define Miniconda installer based on architecture
    if ($ARCH -eq "arm64") {
        $INSTALLER = "Miniconda3-latest-Windows-arm64.exe"
    } else {
        $INSTALLER = "Miniconda3-latest-Windows-x86_64.exe"
    }

    # Download Miniconda Installer
    Invoke-WebRequest -Uri "https://repo.anaconda.com/miniconda/$INSTALLER" -OutFile "$env:USERPROFILE\$INSTALLER"

    # Run the Miniconda installer
    Start-Process -FilePath "$env:USERPROFILE\$INSTALLER" -ArgumentList "/S", "/D=$env:USERPROFILE\miniconda" -Wait

    # Clean up the installer
    Remove-Item "$env:USERPROFILE\$INSTALLER"

    # Add Miniconda to PATH
    $env:PATH = "$env:USERPROFILE\miniconda\Scripts;$env:USERPROFILE\miniconda;$env:PATH"
    Write-Host "Miniconda installed at $env:USERPROFILE\miniconda"
}

# Check if Conda is installed
if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
    Install-Miniconda
} else {
    Write-Host "Conda found: $(Get-Command conda)"
}

# Initialize Conda (in PowerShell)
conda init

# Create or Update Conda environment
$envList = conda env list
if ($envList -match $ENV_NAME) {
    Write-Host "Environment '$ENV_NAME' already exists. Updating environment..."
    conda env update -f $ENV_FILE --name $ENV_NAME
} else {
    Write-Host "Creating new Conda environment from $ENV_FILE..."
    conda env create -f $ENV_FILE
}

# Activate Conda environment
Write-Host "Activating Conda environment..."
conda activate $ENV_NAME

# Install R packages
Write-Host "Running R package installation..."
Rscript scripts/install_packages.R
