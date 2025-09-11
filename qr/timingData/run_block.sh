m=1000
n=${m}
k=512
num_repeat=10
for (( nb=32; nb<=${k}; nb+=32 ))
do
  echo "m=${m} n=${n} k=${k} nb=${nb}"
  echo "Testing file ran ${num_repeat} times"
  for (( l=1; l<=${num_repeat}; l+=1 ))
  do
    ./main.exe -m $m -n $n -k $k -nb $nb
  done
done
