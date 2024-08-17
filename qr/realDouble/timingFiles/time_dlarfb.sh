#!/bin/env bash

# Fix n and k
# n is a decently sized number of columns, and k is the blocksize. We are sticking with
# 32 in this as we are concerned with performance inside dorgqr for dlarfb
# We will later investigate what the performance looks like for larger sizes of nb
for (( m=50000; m<= 200000; m+=2000 ))
do 
    echo "m=$m"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        # n is this small so that we can run it with our innefficient timing method.
        ./timeDlarfb.exe -m $m -n 1000 -k 32 -t
    done
done
