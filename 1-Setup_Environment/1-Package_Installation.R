



####---- Single Cell Pipeline Package Installation ----####


## the following packages had to be installed via github
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remotes::install_github("mojaveazure/seurat-disk")


## Install and load CRAN packages
packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder",
              "Matrix","fields","SeuratDisk","openxlsx")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


## Install and load Bioconductor packages
bioCpacks <- c("SingleR","celldex","SingleCellExperiment","rhdf5")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))
