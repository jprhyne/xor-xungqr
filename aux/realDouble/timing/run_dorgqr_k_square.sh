#!/bin/env bash
dateStr=$(date +"%m_%d_%y")
saveFile="square_orgqr_nb_k${dateStr}.txt"
rm ${saveFile}

m=1000
n=$m

kMin=100

kMax=1000

kDelta=10

numRepeat=5

for (( k=kMin ; k<=kMax ; k+=kDelta )); do
  echo "m=$m:n=$n:k=$k:nb=64" >> ${saveFile}
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n $k | ./time_dorgqr_k.exe >> ${saveFile}
  done
done
