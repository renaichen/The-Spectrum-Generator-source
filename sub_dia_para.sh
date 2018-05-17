#! bin/bash

xleft=()
# for i in `seq 1 3 4`
# do
#     xleft+=($i)
# done
# for i in `seq 0.6 0.1 0.9`
# do
#     xleft+=($i)
# done
# xleft=(0.1 1.0)
xleft=1.0

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

