#!/bin/env bash

# Fix m and vary n=k. 
m=100000
for (( k=1; k<=2048; k*=2 ))
do 
    n=$k
    echo "n=$n"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./overallTest_optBlas.exe -t -m $m -n $n -k $k
    done
done
