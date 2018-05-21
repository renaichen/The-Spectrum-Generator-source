#! /bin/bash
#$ -N Renai-single-multi
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 ~/THE-Spectrum-Generator-Class/src/single_impurity_traj.py $1
