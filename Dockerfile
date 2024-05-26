FROM rocker/r-base:latest


RUN echo "\noptions(shiny.port=3838, shiny.host='0.0.0.0')" >> /usr/lib/R/etc/Rprofile.site

# system libraries of general use
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libarchive-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \ 
    libjpeg-dev \
    libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*




# basic shiny functionality
RUN R -q -e "options(warn=2); install.packages(c('shiny'))"

# install dependencies of the euler app
RUN R -q -e "options(warn=2); install.packages(c('devtools','shinymanager','dplyr','shinyjs','bslib','shinylive','shinycssloaders','DT','mathjaxr','colorspace','ggplot2','gridExtra','cowplot','ggrepel','ggpp','ggh4x','statmod','circlize','tidyverse','tibble','BiocManager'))"
RUN R -q -e "options(warn=2); library(BiocManager)"
RUN R -q -e "options(warn=2); BiocManager::install(version = '3.19')"
RUN R -q -e "options(warn=2); BiocManager::install(c('ComplexHeatmap','standR','limma','scater'))"

# install R code
COPY . /app
WORKDIR /app

EXPOSE 3838

CMD ["R", "-q", "-e", "shiny::runApp('/app')"]