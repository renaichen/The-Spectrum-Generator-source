#! bin/bash

xleft=()
# for i in `seq 1 3 4`
# do
#     xleft+=($i)
# done
for i in `seq 0.7 0.1 1.9`
do
    xleft+=($i)
done
# xleft=1.2

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    sbatch ./sbatch_dia.sh ${xleft[$i-1]} 
   # sleep 18m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done

