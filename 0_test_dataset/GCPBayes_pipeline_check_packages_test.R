# ================================================================================
# 					CHECKING AND INSTALLING A LIST OF PACKAGES
# ================================================================================
# Summary: R script to check a list of packages and install them if they are not available
# It is also print a warning message if any of the packages could not be installed.
# ================================================================================
# Written first by: Yazdan Asgari
# Modified by: Yazdan Asgari
# Initial Creation Date: 02/2023
# Edited Date: 03/2023
# ================================================================================

# Define the list of packages to check and install
packages <- c("vroom", "tidyverse", "data.table", "optparse",
              "devtools", "MASS", "dplyr", "tidyr",
              "tictoc", "GCPBayes",
              "shiny", "datasets", "ggplot2", "gridExtra", "BioCircos", "plotly",
              "genetics.binaRies")

# Create an empty vector to store the names of any packages that could not be installed
missing_packages <- c()

# Loop through the list of packages and install them if they are not available
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)

    # If the package is still not available after installation, add its name to the missing_packages vector
    if (!require(package, character.only = TRUE)) {
      missing_packages <- c(missing_packages, package)
    }
  }
}

# If any packages could not be installed, print a warning message
if (length(missing_packages) > 0) {
  message("Warning: The following packages could not be installed: ", paste(missing_packages, collapse = ", "))
} else {
  message("All the packages are available in your system. Enjoy using GCPBayes Pipeline !")
}

