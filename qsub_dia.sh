#! /bin/bash
#$ -N Renai-diatomic-multi
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 ~/THE-Spectrum-Generator-Class/src/diatomic_white_parallel.py $1
