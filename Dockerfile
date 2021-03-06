## Emacs, make this -*- mode: sh; -*-
## Adapted from rocker/r-base

FROM r-base

LABEL org.label-schema.license="GPLv3.0" \
      maintainer="Malindrie Dharmaratne, <malindrie@gmail.com>" \
      version="0.1.0"

# system libraries of general use
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      libhdf5-dev \
      libv8-dev \
      libssl-dev \
      libcurl4-openssl-dev \
      libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN mkdir -p $HOME/.R/ \
    && echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC\n" >> $HOME/.R/Makevars


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


