#!/bin/env bash
#This file will time driver using test.exe as well as test_optBlas.exe for the reference and optimized BLAS/LAPACK routines resp.

# First, we run make clean to clear out the directory
#make clean

# Next, we run make to compile our 
#make


# Now, we will run test.exe for varying m,n,k files only grabbing the timing
echo "Tall and skinny matrices"
for (( m=15000; m<=155000; m+=10000 ))
do
    echo "m=$m n=32 k=32"
    echo "Testing file ran 5 times"
    for (( l=1; l<=5; l+=1 ))
    do
        ./overallTest_optBlas.exe -m $m -n 32 -k 32 -t
    done
done
