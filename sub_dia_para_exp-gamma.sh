#! bin/bash

xleft=()
# for i in `seq 1 1 5`
# do
#     xleft+=($i)
# done

for j in `seq 0.65 0.025 1.0`
do
    xleft+=($j)
done

xleft+=(2 3 4 5)
# xleft+=(0.8 0.5)
# xleft=0.8

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    qsub ~/THE-Spectrum-Generator-Class/src/qsub_dia_exp-gamma.sh ${xleft[$i-1]} 
    # qsub qsub_single.sh ${xleft[$i-1]} 
   # sleep 18m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done
