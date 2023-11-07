# Use the rocker/r-ver image for R 4.3.1
FROM rocker/r-ver:4.3.1

# Maintainer information
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"

# Update and install required system packages and clean up
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev \
    libudunits2-dev \
    libxml2-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk40 && \
    rm -rf /var/lib/apt/lists/*

# Set library paths for R packages
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib:/usr/local/lib/R/lib

# Install necessary R packages
RUN R -e 'install.packages(c("shiny", "shinythemes", "shinyjqui", "pheatmap", "RColorBrewer", "umap", "shinyjs", "cowplot", "patchwork", "ggdendro", "factoextra", "dplyr", "DT", "viridis", "readr", "ggrepel", "ggVennDiagram", "ggtree", "shinycssloaders", "stringr", "tools", "plotly", "reshape2", "ggpubr", "gridExtra", "scales", "tidyverse"))'

# Install Bioconductor packages
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' \
  && R -e 'BiocManager::install(c("clusterProfiler", "GSVA", "slingshot", "SingleCellExperiment"), ask = FALSE)'

# Install other R packages
RUN R -e 'install.packages(c("remotes", "juicyjuice", "xml2", "gt", "gtsummary", "gridtext", "ggtext", "survminer"))'

# Set the Shiny port and host
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" 

# Create a user for running the app
RUN addgroup --system app && adduser --system --ingroup app app

# Copy the R Shiny app code to the container
COPY ./app/ /app

# Set user permissions
RUN chown app:app -R /app

# Expose the port
EXPOSE 3838

# Start the Shiny app
USER app
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
