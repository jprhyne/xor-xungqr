#!/bin/env bash

# Fix m
m=100000

for (( n=100; n<=1000; n+=50 ))
do
    echo "n=$n"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeDorgqr.exe -m $m -n $n -k $n -t
    done
done
