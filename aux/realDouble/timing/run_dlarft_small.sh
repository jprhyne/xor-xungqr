#!/bin/env bash
rm small.out

m=10000
nMin=10

nMax=1000

nDelta=10

numRepeat=5

for (( n=nMin ; n<=nMax ; n+=nDelta )); do
  echo "m=$m:n=$n" >> small.out
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n | ./time_dlarft.exe >> small.out
  done
done
