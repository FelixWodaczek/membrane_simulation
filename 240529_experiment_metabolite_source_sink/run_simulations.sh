#!/bin/bash -l

#SBATCH --job-name=Cntrl6
#SBATCH --output="output/output-%A_%a_%x.out"
#SBATCH --exclude=zeta[243-263]
#SBATCH --nodes 1
#SBATCH -c 4
#SBATCH --time=240:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=maitane.munoz-basagoiti@ist.ac.at
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --export=None

env > 'env_file.dat'
unset SLURM_EXPORT_ENV

#Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#load the respective software module you intend to use
# here one would need to load the sconda environment for Trilmp, here called Trimenv
module purge
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
module load openmpi

conda activate Trienv

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}

path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bashdirs/directories_7.dat)
cd ${path}

#run the respective binary through SLURM's srun
srun --cpu_bind=verbose python3 launch.py