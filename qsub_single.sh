#! /bin/bash
#$ -N Renai-single-multi
#$ -V
#$ -cwd
#$ -pe openmp 46

python2 ~/THE-Spectrum-Generator-Class/src/single_white_harmonic_NoTrajHold.py $1
