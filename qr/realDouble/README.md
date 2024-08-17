# orgqr

This repository will hold my modifications on the dorgqr.f file from release 3.11.0.

# Compilation instructions

1. Rename make.inc.example to make.inc
2. Fill in the paths to LAPACK libraries 
3. Run `make`

# Testing instructions

1. Follow the compilation instructions
2. Run `./driver > results.txt`
3. Run `./checkMax.py`
4. Observe results for maximum relative representation and orthogonality errors

# Files

* checkMax.py 
    * This file is a short python script that helps parse the output of the driver file
* dorgqr_mod.f
    * File used to test changing the NB parameter in the original code structure
* driver.c
    * Quick test driver file to run test.exe with many parameters for a more exhaustive test setup
* LICENSE
    * License for this repository
* Makefile
    * Makefile to compile executables in this repository
* make.inc.example
    * Example make.inc file to give the locations of lapack and blas libraries
* my_dlarfb.f
    * Modified dlarfb for the first iteration used in the dorgqr function. Currently not used
    * TODO: Remove this file or properly call it inside my_dorgqr.f
* my_dorgqr.f
    * Modified dorgqr
* my_dorgqr_v1.f
    * Modified dorgqr that only takes into account of the initial Identity block
* README.md
    * This file
* test.c
    * Test file that tests the representation and orthogonality of the output from my_dorgqr.f

