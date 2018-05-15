#! bin/bash

xleft=()
for i in `seq 1 3 4`
do
    xleft+=($i)
done

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    echo $i >> index.txt
    # qsub qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_single.sh ${xleft[$i-1]} 
   # sleep 5m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done

