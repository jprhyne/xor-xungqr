#!/bin/env bash

# Fix m and vary n=k
m=20000
for (( n=5000; n<=20000; n+=2000 ))
do 
    k=$n
    echo "n=$n"
    # Since n = k at this point, we are only testing the case of
    # computing Q1 and Q2 simultaneously
    echo "Computing Q1 and Q2 simultaneously"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeAll.exe -t -m $m -n $n -k $k -v 1
    done
done
