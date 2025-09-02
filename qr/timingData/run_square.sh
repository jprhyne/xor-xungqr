#!/bin/env bash

# Fix m and n, vary k
m=10000
n=${m}
num_repeat=10
for (( k=1000; k<=${m}; k+=1000 ))
do
  echo "m=${m} n=${n} k=${k}"
  echo "Testing file ran ${num_repeat} times"
  for (( l=1; l<=${num_repeat}; l+=1 ))
  do
    ./main.exe -m $m -n $n -k $k
  done
done
