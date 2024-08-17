#!/bin/env bash
#This file will time driver using test.exe as well as test_optBlas.exe for the reference and optimized BLAS/LAPACK routines resp.

# First, we run make clean to clear out the directory
#make clean

# Next, we run make to compile our 
#make


# Now, we will run test.exe for varying m,n,k files only grabbing the timing
echo "Tall and skinny matrices"
for (( m=500000; m<=500000; m+=200000 ))
do
    for (( n=10; n<=100; n+=10 ))
    do
        echo "m=$m n=$n k=$n"
        echo "Testing file ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./overallTest_optBlas.exe -m $m -n $n -k $n -t
        done
    done
done

echo "T-Shorter matrices"
for (( m=35000; m<=335000; m+= 100000 ))
do
    for (( n=32; n<=532; n+=100 ))
    do
        echo "m=$m n=$n k=$n"
        echo "Testing file ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./overallTest_optBlas.exe -m $m -n $n -k $n -t
        done
    done
done

