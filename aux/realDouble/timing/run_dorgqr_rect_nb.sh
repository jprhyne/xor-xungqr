#!/bin/env bash
dateStr=$(date +"%m_%d_%y")
saveFile="rect_orgqr_nb_${dateStr}.txt"
rm ${saveFile}

m=10000
n=1000
k=1000

nbMin=10
nbMax=1000
nbDelta=10

numRepeat=5

for (( nb=nbMin ; nb<=nbMax ; nb+=nbDelta )); do
  echo "m=$m:n=$n:k=$k:nb=${nb}" >> ${saveFile}
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n $k ${nb} | ./time_dorgqr.exe >> ${saveFile}
  done
done
