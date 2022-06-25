#!/bin/bash
# Script for installing stuff on CentOS
# Requires sudo privileges

#########
# Useful stuff
function python3 {
    yum -y install yum-utils
    yum-builddep python
    curl -O https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tgz
    tar xf Python-3.5.0.tgz
    cd Python-3.5.0
    ./configure
    make
    make install
}

# Pip and get some packages for Python
function pypackages {
    yum -y install ipython
    yum -y install python-pip
    pip install --upgrade pip

    pip install pandas
    pip3 install pandas
    pip install astropy
    pip3 install astropy
    pip install seaborn
    pip3 install seaborn
}

#########
# Prerequisites for XSPEC
# Install some stuff for XSPEC
# Set some variables for XSPEC
function prep_xspec {
    yum -y install ncurses-devel
    yum -y install libXt-devel
    yum -y install gcc-c++
    yum -y install gcc gcc-gfortran
    yum -y install perl-ExtUtils-MakeMaker
    yum -y install python-devel
    yum -y install libpng-devel
    
    export CC=/usr/bin/gcc
    export CXX=/usr/bin/g++
    export FC=gfortran
    export PERL=perl
    export PYTHON=/usr/bin/python
}

# Install XSPEC, requires the folder with the installation failes
function setup_xspec {
    gunzip -c heasoft-6.19.src.tar.gz | tar xf -
    cd heasoft-6.19/BUILD_DIR
    ./configure
    make
    make install
}
