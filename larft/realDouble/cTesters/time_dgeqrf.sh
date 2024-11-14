minM=5000
maxM=50000
incM=1000
echo "Tall current case of 32 block size"
for (( m=$minM; m<$maxM; m+=$incM ))
do
    echo "m=$m n=$minM nb=32"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDgeqrf.exe -m $m -n $minM -nb 32
    done
done
echo "Wide current case of 32 block size"
minN=10
maxN=2000
incN=100
for (( n=$minN; n<$maxN; n+=$incN ))
do
    echo "m=$minM n=$n nb=32"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDgeqrf.exe -m $minM -n $n -nb 32
    done
done
echo "Trying new block sizes"
# fix m and vary n
m=10000
maxN=8192
for (( n=16; n<=maxN ; n*=2 ))
do
    echo "m=$m n=$m nb=$n"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDgeqrf.exe -m $m -n $m -nb $n
    done
done
