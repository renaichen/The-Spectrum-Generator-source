#! bin/bash

xleft=()
# for i in `seq 1 1 5`
# do
#     xleft+=($i)
# done

# for j in `seq 0.1 0.1 0.9`
# do
#     xleft+=($j)
# done

xleft+=(2 3)
# xleft=4

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    qsub ~/THE-Spectrum-Generator-Class/src/qsub_dia.sh ${xleft[$i-1]} 
    # qsub qsub_single.sh ${xleft[$i-1]} 
   # sleep 18m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done

