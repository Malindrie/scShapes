Bootstrap: docker
From: ubuntu:16.04

%help
  This container runs R.

%labels
  Maintainer Malindrie Dharmaratne.

%apprun R
  exec R "${@}"

%apprun Rscript
  exec Rscript "${@}"

%runscript
  exec R "${@}"

%post
    mkdir -p /scratch/global /scratch/local /rcc/stor1/refdata /rcc/stor1/projects /rcc/stor1/depts /extR/library1 /extR/library2
    apt-get update
    apt-get -y install \
        wget \
        build-essential \
        software-properties-common \
        apt-transport-https \
        locales \
        libv8-dev \
        libxml2-dev \
        libhdf5-serial-dev
    apt-get -y install libcurl4-openssl-dev libssl-dev
    echo "LC_ALL=en_US.UTF-8" >> /etc/environment
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
    echo "LANG=en_US.UTF-8" > /etc/locale.conf
    locale-gen en_US.UTF-8
    echo 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' >> /etc/apt/sources.list
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 51716619E084DAB9
    add-apt-repository -y "ppa:marutter/rrutter3.5"
    add-apt-repository -y "ppa:marutter/c2d4u3.5"
    apt-get update
    apt-get -y install --no-install-recommends --allow-unauthenticated\
        r-base \
        r-base-core \
        r-base-dev
    apt-get clean

  mkdir -p $HOME/.R/
    echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC\n" >> $HOME/.R/Makevars

  Rscript -e "install.packages('Matrix')"
  Rscript -e "install.packages('stats')"
  Rscript -e "install.packages('parallelly')"
  Rscript -e "install.packages('pscl')"
  Rscript -e "install.packages('emdbook')"
  Rscript -e "install.packages('future')"
  Rscript -e "install.packages('future.apply')"
  Rscript -e "install.packages('pscl')"
  Rscript -e "install.packages('VGAM')"
  Rscript -e "install.packages('MASS')"
  Rscript -e "install.packages('magrittr')"
  Rscript -e "install.packages('utils')"
  Rscript -e "install.packages('dgof')"
  Rscript -e "install.packages('devtools')"
  Rscript -e "devtools::install_github('Malindrie/scShapes', dep=FALSE, build_vignettes=TRUE)"

  rm -rf /var/lib/apt/lists/*

  echo "options(repos = c(CRAN = 'https://cloud.r-project.org/'), download.file.method = 'libcurl')" >> /usr/lib/R/etc/Rprofile.site
