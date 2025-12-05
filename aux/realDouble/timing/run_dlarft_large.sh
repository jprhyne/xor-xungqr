#!/bin/env bash

# Large test cases
m=10000
nMin=1000

nMax=10000

nDelta=500

numRepeat=5

rm large.out
for (( n=nMin ; n<=nMax ; n+=nDelta )); do
  echo "m=$m:n=$n" >> large.out
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n | ./time_dlarft.exe >> large.out
  done
done


