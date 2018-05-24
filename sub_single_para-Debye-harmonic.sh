#! bin/bash

xleft=()
# for i in `seq 1 3 4`
# do
#     xleft+=($i)
# done
xleft=24
# xleft=(500 2500 5000)

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    qsub ~/THE-Spectrum-Generator-Class/src/qsub_single-Debye-harmonic.sh ${xleft[$i-1]} 
   # sleep 5m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done

