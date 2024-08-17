#!/bin/env bash

# Fix m and n vary k. 
m=10000
n=$m
for (( k=1000; k<=$m; k+=200 ))
do 
    echo "k=$k"
    echo "Computing Q1 and Q2 simultaneously"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./test.exe -v 8 -t -m $m -n $n -k $k
    done
    echo "Computing Q1 and Q2 separately"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./test.exe -v 9 -t -m $m -n $n -k $k
    done
done
