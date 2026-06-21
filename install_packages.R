# install_packages.R

# 1. Define standard CRAN packages
cran_packages <- c("tidyverse", "ape", "phytools", "progressr", "future", "furrr") 
install.packages(cran_packages, repos = "https://r-project.org")

# 2. Install 'remotes' to handle GitHub installations
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://r-project.org")
}

# 3. Define and install GitHub packages (Format: "developer/repository")
github_packages <- c(
  "GabrielNakamura/Herodotools",
  "vanderleidebastiani/daee",
  "nmatzke/BioGeoBEARS"
)
remotes::install_github(github_packages, upgrade = "never")
