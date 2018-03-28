#!/bin/bash
#
# setup and configure the HapCHAT core phasing algorithm
#
base="$(git rev-parse --show-toplevel)"

hapchatcore=$base/build/hapchat_core

#software/hapchat/hapchat_core

# The HapCHAT core phasing algorithm
#----------------------------------------------------------------------
echo "setting up the HapCHAT core phasing algorithm .."
echo

# add a build directory and build the source
function build_hapchatcore {
    cd $base
    echo "building source with cmake/make .."
    mkdir -p build
    cd build
    cmake ../src
    make -j 8
}

if [ -d build ] ; then
    echo "deleting old instance .."
    rm -rf build
fi
    
build_hapchatcore

echo "setup of the HapCHAT core phasing algorithm is complete."

# finish
#----------------------------------------------------------------------
echo
echo "The HapCHAT core phasing algorithm can now be executed as:"
echo
echo "   $hapchatcore"
echo
