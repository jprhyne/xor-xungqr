#!/bin/env bash
dateStr=$(date +"%m_%d_%y")
saveFile="small_larft_${dateStr}.txt"
rm ${saveFile}


m=10000
nMin=10

nMax=1000

nDelta=10

numRepeat=5

for (( n=nMin ; n<=nMax ; n+=nDelta )); do
  echo "m=$m:n=$n" >> ${saveFile}
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n | ./time_dlarft.exe >> ${saveFile}
  done
done
