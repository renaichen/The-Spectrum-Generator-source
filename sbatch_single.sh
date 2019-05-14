#! /bin/bash
#SBATCH --job-name="Renai-single-mult"
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-node=24
#SBATCH -A sds154

module load python
module load scipy

# python ~/DebyeBath/The-Spectrum-Generator-source/single_white_parallel.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/single_impurity_singletraj.py $1
python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/single_impurity_traj.py $1
