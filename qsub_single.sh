#! /bin/bash
#$ -N Renai-multi-single-test
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 single_impurity_singletraj.py $1
