#!/bin/env bash
# Fix m and vary n=k
m=20000
for (( n=5000; n<=20000; n+=2000 ))
do 
    k=$n
    echo "m=$m, n=$n"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeDorgqrVsDorglq.exe -m $m -n $n -k $k
    done
done
