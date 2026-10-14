# R packages for notebooks 1 and 3 (original analysis: R 4.3.x, Jupyter IRkernel).
# Record exact versions after running the notebooks with: sessionInfo()

cran_pkgs <- c("dplyr", "openxlsx", "data.table", "stringr",
               "Seurat", "HGNChelper", "IRkernel")
install.packages(setdiff(cran_pkgs, rownames(installed.packages())))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("flowCore", update = FALSE)

# Register the R kernel for Jupyter
IRkernel::installspec(name = "ir", displayname = "R")

# ScType functions are sourced at runtime in notebook 3 from
# https://github.com/IanevskiAleksandr/sc-type (R/gene_sets_prepare.R, R/sctype_score_.R)
