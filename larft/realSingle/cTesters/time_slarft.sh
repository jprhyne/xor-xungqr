for (( m=500000; m<5000000; m+=200000 ))
do
    echo "m=$m"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeSlarft.exe -t -m $m
    done
done

# fix m and vary n
m=500000
for (( n=1; n<=2048 ; n*=2 ))
do
    echo "n=$n"
    echo "Teting file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeSlarft.exe -t -m $m -n $n
    done
done
