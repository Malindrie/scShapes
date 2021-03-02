BootStrap: docker
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
  apt-get update
  apt-get install -y apt-transport-https apt-utils software-properties-common
  apt-get install -y add-apt-key
  export DEBIAN_FRONTEND=noninteractive
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  apt-get install -y tzdata
  dpkg-reconfigure --frontend noninteractive tzdata

  #add CRAN/Ubuntu repo, add key, then refresh
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
  apt-get update

  apt-get install -y wget nano
  apt-get install -y libblas3 libblas-dev liblapack-dev liblapack3 curl
  apt-get install -y gcc fort77 aptitude
  aptitude install -y g++
  aptitude install -y xorg-dev
  aptitude install -y libreadline-dev
  aptitude install -y gfortran
  gfortran --version
  apt-get install -y libssl-dev libxml2-dev libpcre3-dev liblzma-dev libbz2-dev libcurl4-openssl-dev
  apt-get install -y libhdf5-dev hdf5-helpers libmariadb-client-lgpl-dev
  apt-get -y install libcurl4-openssl-dev libssl-dev
  apt-get update
  apt-get -y install --no-install-recommends --allow-unauthenticated\
    r-base \
    r-base-core \
    r-base-dev \
    r-recommended \
    r-base-html \
    r-doc-html \
    r-cran-devtools
  apt-get clean

  mkdir -p $HOME/.R/
    echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC\n" >> $HOME/.R/Makevars

  Rscript -e "install.packages('parallelly')"
  Rscript -e "install.packages('pscl')"
  Rscript -e "install.packages('emdbook')"
  Rscript -e "install.packages('future')"
  Rscript -e "install.packages('future.apply')"
  Rscript -e "install.packages('pscl')"
  Rscript -e "install.packages('dgof')"
  Rscript -e "devtools::install_github('Malindrie/scShapes', dep=FALSE, build_vignettes=TRUE)"

  rm -rf /var/lib/apt/lists/*

  echo "options(repos = c(CRAN = 'https://cloud.r-project.org/'), download.file.method = 'libcurl')" >> /usr/lib/R/etc/Rprofile.site
