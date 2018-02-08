#!/bin/bash
#
# setup and configure WhatsHap and HapCHAT
#
base="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd )"

whatshap=$base/software/whatshap/venv/bin/whatshap
hapchat=$base/software/hapchat/hapchat_core

# WhatsHap
#----------------------------------------------------------------------
echo "setting up WhatsHap .."

# clone git repo
function clone_whatshap {
    cd $base
    echo "cloning WhatsHap git repository .."
    mkdir -p software
    cd software
    git clone --branch HapCHAT https://bitbucket.org/whatshap/whatshap.git
}

# setup virtual environment
function venv_whatshap {
    cd $base
    echo "setting up virtual environment .."
    cd software/whatshap
    virtualenv -p python3 venv
    venv/bin/pip3 install Cython nose tox
    venv/bin/pip3 install -e .
}

if [ -e $whatshap ] ; then
    echo "WhatsHap excecutable already exists"
else
    clone_whatshap
    venv_whatshap
fi

echo "setup of WhatsHap is complete."

# HapCHAT
#----------------------------------------------------------------------
echo "setting up the HapCHAT core phasing algorithm .."

# add a build directory and build the source
function build_hapchat {
    cd $base
    echo "building source with cmake/make .."
    mkdir -p software
    cd software
    mkdir -p hapchat
    cd hapchat
    cmake ../../src
    make -j 8
}

if [ -e $hapchat ] ; then
    echo "HapCHAT executable already exists"
else
    build_hapchat
fi

echo "setup of HapCHAT core phasing algorithm is complete."

# finish
#----------------------------------------------------------------------

echo
echo "setup of both tools has been successful." 
echo
echo "WhatsHap can now be executed as:"
echo
echo "   $whatshap"
echo
echo "HapCHAT can now be executed as:"
echo
echo "   $hapchat"
echo
