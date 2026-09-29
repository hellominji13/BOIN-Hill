# hash:sha256:e36078db4bad8cd6a8690dbffe0107cb7a7d39693a552cb0f7700d5b7e17773a
FROM registry.codeocean.com/codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        jags=4.3.2-1.2204.0 \
        libgsl-dev=2.7.1+dfsg-3 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("rjags", "4-17")' \
    && Rscript -e 'remotes::install_version("MASS", "7.3-66")' \
    && Rscript -e 'remotes::install_version("shiny", "1.14.0")' \
    && Rscript -e 'remotes::install_version("rhandsontable", "0.3.8")' \
    && Rscript -e 'remotes::install_version("LaplacesDemon", "16.1.8")' \
    && Rscript -e 'remotes::install_version("ggplot2", "4.0.3")' \
    && Rscript -e 'remotes::install_version("readxl", "1.5.0.1")' \
    && Rscript -e 'remotes::install_version("future", "1.76.0")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.17")' \
    && Rscript -e 'remotes::install_version("foreach", "1.5.2")' \
    && Rscript -e 'remotes::install_version("BayesLogit", "2.4")' \
    && Rscript -e 'remotes::install_version("VGAM", "1.1-14")' \
    && Rscript -e 'remotes::install_version("BOIN", "2.7.2")' \
    && Rscript -e 'remotes::install_version("UBCRM", "1.0.3")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.5")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.2.1")' \
    && Rscript -e 'remotes::install_version("readr", "2.2.0")' \
    && Rscript -e 'remotes::install_version("copula", "1.1-7")' \
    && Rscript -e 'remotes::install_version("gsl", "2.1-9")'

COPY postInstall /
RUN /postInstall
