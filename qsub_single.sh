#! /bin/bash
#$ -N Renai-multi-single-test
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 ~/THE-Spectrum-Generator-Class/src/single_white_parallel.py $1
