#!/bin/env bash
# This file is how I am compiling all subdirectories. In order to use this file, you must first do either
#   1) copy this file into each sub folder (ie ./realDouble/test/)
#   2) in each sub folder (ie ./realDouble/test) make a symbolic link to this file by doing `ln -s ../../compile.sh`
# Remove any generated archive that may be present
rm *.a
# Compile the lapack source files
make -C ${GIT_REPO_LOC}/myLapackFork/SRC -j8
# Add my functions to my library
myLapackLib
# clean up this directory
make clean
# compile here
make
