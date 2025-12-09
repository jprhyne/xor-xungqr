#!/bin/env bash
dateStr=$(date +"%m_%d_%y")
saveFile="square_larfb_${dateStr}.txt"
rm ${saveFile}

m=10000
n=$m

kMin=10

kMax=1000

kDelta=10

numRepeat=5

for (( k=kMin ; k<=kMax ; k+=kDelta )); do
  echo "m=$m:n=$n:k=$k" >> ${saveFile}
  for (( i=1 ; i<=numRepeat ; i++ )); do
    echo $m $n $k | ./time_dlarfb0c2.exe >> ${saveFile}
  done
done
