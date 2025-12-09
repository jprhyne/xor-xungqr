#!/bin/env bash
dateStr=$(date +"%m_%d_%y")
saveFile="rect_orgqr_${dateStr}.txt"
rm ${saveFile}

m=10000

kMin=100

kMax=10000

kDelta=100

numRepeat=5

for (( k=kMin ; k<=kMax ; k+=kDelta )); do
  n=$k
  echo "m=$m:n=$n:k=$k:nb=64" >> ${saveFile}
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n $k 64 | ./time_dorgqr.exe >> ${saveFile}
  done
done
