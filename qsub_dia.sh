#! /bin/bash
#$ -N Renai-multi-diatomic-test
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 ~/THE-Spectrum-Generator-Class/src/diatomic_impurity_traj.py $1
