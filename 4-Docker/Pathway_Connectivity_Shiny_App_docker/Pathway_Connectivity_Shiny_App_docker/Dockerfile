# Use the rocker/r-ver image for R 4.3.1
FROM rocker/r-ver:4.3.1

# Maintainer information
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    cmake \
    libxml2-dev \
    libglpk40 \
    && rm -rf /var/lib/apt/lists/*

# Set library paths for R packages
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib:/usr/local/lib/R/lib

# Install R packages non-interactively to the user-level library
RUN R -e 'install.packages(c("shiny", "shinyjqui", "pheatmap", "RColorBrewer", "ggdendro", "factoextra", "dplyr", "DT", "viridis", "readr", "ggrepel", "shinycssloaders", "stringr", "tools", "plotly", "reshape2", "ggpubr", "gridExtra"))'

# Install Bioconductor packages to the user-level library
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' \
  && R -e 'BiocManager::install("clusterProfiler", ask = FALSE)'

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
