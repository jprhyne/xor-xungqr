#!/bin/env bash

# Fix m and n, vary k
m=1000
n=${m}
num_repeat=10
for (( k=200; k<=${m}; k+=100 ))
do
  nb=${k}
  echo "m=${m} n=${n} k=${k} nb=${nb}"
  echo "Testing file ran ${num_repeat} times"
  for (( l=1; l<=${num_repeat}; l+=1 ))
  do
    ./main.exe -m $m -n $n -k $k
  done
done
