#! /bin/bash
#SBATCH --job-name="Renai-diatomic-mult"
#SBATCH --time=02:00:00
#SBATCH --ntasks-per-node=24
#SBATCH -A sds154

module load python
module load scipy

# python ~/DebyeBath/The-Spectrum-Generator-source/diatomic_white_parallel.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/diatomic_impurity_traj.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_impurity_traj.py $1
python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_impurity_traj_harmonic.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_white_harmonic_NoTrajHold.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_white_harmonic_NoTrajHold-gamma.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_white_exponential_NoTrajHold.py $1
# python ~/DebyeBath/The-Spectrum-Generator-source/NoTrajHold/diatomic_white_exponential_NoTrajHold-gamma.py $1
