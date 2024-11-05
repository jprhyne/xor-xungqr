#minM=500000
#maxM=5000000
#incM=200000
minM=5000
maxM=50000
incM=2000
for (( m=$minM; m<$maxM; m+=$incM ))
do
    echo "m=$m"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDgeqrf.exe -m $m -n $m
    done
done

# fix m and vary n
m=$minM
maxN=2048
for (( n=1; n<=maxN ; n*=2 ))
do
    echo "n=$n"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDlarft.exe -t -m $m -n $n
    done
done
