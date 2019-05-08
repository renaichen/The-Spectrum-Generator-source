#! /bin/bash
#SBATCH --job-name="Renai-diatomic-mult"
#SBATCH --time=00:10:00
#SBATCH --ntasks-per-node=24
#SBATCH -A sds154

module load python
module load scipy

python ~/DebyeBath/The-Spectrum-Generator-source/diatomic_white_parallel.py $1
