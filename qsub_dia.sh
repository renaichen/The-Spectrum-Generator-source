#! /bin/bash
#$ -N Renai-multi-diatomic-test
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 diatomic_impurity_traj.py $1
