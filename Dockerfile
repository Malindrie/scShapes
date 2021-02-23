## Emacs, make this -*- mode: sh; -*-
## Adapted from rocker/r-base

FROM r-base

LABEL org.label-schema.license="GPLv3.0" \
      maintainer="Malindrie Dharmaratne, <malindrie@gmail.com>" \
      version="0.1.0"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev


# install R packages required
RUN R --slave -e 'install.packages("Matrix")'
RUN R --slave -e 'install.packages("stats")'
RUN R --slave -e 'install.packages("parallelly")'
RUN R --slave -e 'install.packages("future")'
RUN R --slave -e 'install.packages("future.apply")'
RUN R --slave -e 'install.packages("pscl")'
RUN R --slave -e 'install.packages("VGAM")'
RUN R --slave -e 'install.packages("dgof")'
RUN R --slave -e 'install.packages("MASS")'
RUN R --slave -e 'install.packages("emdbook")'
RUN R --slave -e 'install.packages("magrittr")'
RUN R --slave -e 'install.packages("utils")'
RUN R --slave -e 'install.packages("devtools")'
RUN R --slave -e 'devtools::install_github("Malindrie/scShapes", dep = FALSE, build_vignettes = TRUE)'


CMD ["/bin/bash"]


