n=5000
minK=10
maxK=5100
incK=100
for (( k=$minK; k<$maxK; k+=$incK ))
do
    echo "n=$n k=$k"
    echo "Testing file ran 10 times"
    for (( l=1; l<= 10; l+=1 ))
    do
        ./timeDlarft.exe -t -m $n -n $k
    done
done
