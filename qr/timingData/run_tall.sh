#!/bin/env bash

# Fix n and k, vary m
n=1000
k=${n}
max_m=10000
num_repeat=10
for (( m=1000; m<=${max_m}; m+=1000 ))
do
  echo "m=${m} n=${n} k=${k}"
  echo "Testing file ran ${num_repeat} times"
  for (( l=1; l<=${num_repeat}; l+=1 ))
  do
    ./main.exe -m $m -n $n -k $k
  done
done
