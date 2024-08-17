#!/bin/env bash
#This file will time driver using test.exe as well as test_optBlas.exe for the reference and optimized BLAS/LAPACK routines resp.

# First, we run make clean to clear out the directory
#make clean

# Next, we run make to compile our 
#make


# Now, we will run test.exe for varying m,n,k files only grabbing the timing
echo "Tall and skinny matrices"
for (( m=200000; m<=1000000; m+=200000 ))
do
    for (( n=50; n<=200; n+=50 ))
    do
        for (( k=50; k<=n; k+=50 ))
        do
            echo "m=$m n=$n k=$k"
            echo "Reference dorgqr ran 30 times"
            for (( l=1; l<=30; l+=1 ))
            do 
                ./timeDorgqr.exe -m $m -n $n -k $k -t
            done
            echo "my_dorgqr ran 30 times"
            for (( l=1; l<=30; l+=1 ))
            do
                ./test_v3.exe -m $m -n $n -k $k -t
            done 
        done
    done
done
echo "Square matrices"
for (( m=1000; m<=5000; m+=1000 ))
do
    n=$m
    for (( k=100; k<=n; k+=100 ))
    do
        echo "m=$m n=$n k=$k"
        echo "Reference dorgqr ran 30 times"
        for (( l=1; l<=30; l+=1 ))
        do 
            ./timeDorgqr.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr ran 30 times"
        for (( l=1; l<=30; l+=1 ))
        do
            ./test_v3.exe -m $m -n $n -k $k -t
        done 
    done
done
