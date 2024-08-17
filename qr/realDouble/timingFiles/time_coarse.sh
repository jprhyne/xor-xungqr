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
        echo "AOCL ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./timeDorgqr.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v0 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v0.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v1 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v1_opt.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v2 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v2.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v3 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v3.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v4 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v4.exe -m $m -n $n -k $n -t
        done
        echo "my_dorgqr_v5 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v5.exe -m $m -n $n -k $n -t
        done
    done
done
echo "Square matrices"
for (( m=5000; m<=5000; m+=1000 ))
do
    n=$m
    for (( k=1000; k<=n; k+=1000 ))
    do
        echo "m=$m n=$n k=$k"
        echo "AOCL ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./timeDorgqr.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v0 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v0.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v1 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v1_opt.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v2 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v2.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v3 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v3.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v4 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v4.exe -m $m -n $n -k $k -t
        done
        echo "my_dorgqr_v5 ran 5 times"
        for (( l=1; l<=5; l+=1 ))
        do
            ./test_v5.exe -m $m -n $n -k $n -t
        done
    done
done
