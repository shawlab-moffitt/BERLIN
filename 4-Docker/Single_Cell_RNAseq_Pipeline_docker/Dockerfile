# Use the rocker/r-ver image for R 4.3.1
FROM rocker/r-ver:4.3.1

# Maintainer information
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"

# Set the working directory inside the container
WORKDIR /app

# Install required system packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2 \
    libxml2-dev \
    libudunits2-dev \
    zlib1g-dev \
    patch 
    
# Install additional dependencies

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk40

RUN apt-get install -y \
    libz-dev \
    libhdf5-dev

# Set library paths for R packages
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib:/usr/local/lib/R/lib


# Install and load the remotes package to install packages from GitHub
RUN R -e "install.packages('remotes')"

# Install devtools
RUN R -e "remotes::install_github('r-lib/devtools')"

# Install hdf5r
RUN R -e "install.packages('hdf5r')"

# Install packages from GitHub
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN R -e "devtools::install_github('satijalab/seurat-data')"

# Install 'celldex' from GitHub
RUN R -e "remotes::install_github('LTLA/celldex')"

# Install Seurat version 4.4.0 and Matrix 1.6.1 from GitHub without updating other packages
RUN R -e "devtools::install_version('Seurat', version = '4.4.0', repos='https://cran.rstudio.com/', dependencies=TRUE, upgrade='never')"
RUN R -e "devtools::install_version('Matrix', version = '1.6.1', repos='https://cran.rstudio.com/')"

# Install R packages from CRAN and Bioconductor
RUN R -e "install.packages(c('dplyr', 'patchwork', 'readr', 'DoubletFinder', 'fields', 'openxlsx'), repos='https://cran.rstudio.com/')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('SingleR', 'celldex', 'SingleCellExperiment', 'rhdf5'), update = TRUE, ask = FALSE)"

# Copy R script and data files into the container
COPY app/Single_Cell_RNAseq_Analysis.R /app/
COPY app/Input/Count_File_Example_GSE116256_AML921A-D0.txt.gz /app/Input/
COPY app/Input/Meta_File_Example_GSE116256_AML921A-D0.txt.gz /app/Input/
COPY app/R/write.h5.R /app/R/

# Run R script when the container starts
CMD ["Rscript", "Single_Cell_RNAseq_Analysis.R"]
