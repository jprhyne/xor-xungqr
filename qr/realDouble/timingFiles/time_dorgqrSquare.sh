#!/bin/env bash

# Fix m
m=3000
n=3000
for (( k=100; k<=$n; k+=100 ))
do
    echo "k=$k"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeDorgqr.exe -m $m -n $n -k $k -t
    done
done
