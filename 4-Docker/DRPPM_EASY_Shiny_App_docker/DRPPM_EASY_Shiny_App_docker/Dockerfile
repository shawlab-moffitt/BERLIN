# Use the rocker/r-ver image for R 4.3.1
FROM rocker/r-ver:4.3.1

# Maintainer information
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"

# Install system dependencies, including libxml2-dev
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    cmake \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk40

# Set library paths for R packages
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib:/usr/local/lib/R/lib

# Install R packages non-interactively to the user-level library
RUN R -e 'install.packages(c("shinythemes", "shinyjqui", "pheatmap", "BiocManager", "RColorBrewer", "dplyr", "DT", "ggplot2", "ggpubr", "reshape2", "tibble", "viridis", "scales", "plotly", "readr", "enrichR", "ggrepel", "tidyr", "tools", "shinycssloaders"))'

# Install Bioconductor packages to the user-level library
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' \
  && R -e 'BiocManager::install(c("GSVA","XML", "annotate", "GSEABase", "clusterProfiler", "limma", "enrichplot"), ask = FALSE)'

# Set the Shiny port and host
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" 

# Create a user for running the app
RUN addgroup --system app \
    && adduser --system --ingroup app app

# Copy your R code and data into the container
COPY ./app/ /app

# Set ownership of app files
RUN chown app:app -R /app

# Expose the port that Shiny app will listen on
EXPOSE 3838

# Define the command to run the Shiny app
USER app
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
