#!/bin/bash
#
#----------------------------------------------------------------
# running a multiple independent jobs
#----------------------------------------------------------------
#
#  Defining options for slurm how to run
#----------------------------------------------------------------
#
#SBATCH --job-name=membraneSimulation
#SBATCH --output=logs/pa_%A-%a.log 
#        %A and %a are placeholders for the jobid and taskid, resp.
#
#Number of CPU cores to use within one node
#SBATCH -c 8
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=60:00:00
#
#Define the amount of RAM used by your job in GigaBytes
#In shared memory applications this is shared among multiple CPUs
#SBATCH --mem=10G
#
#Do not requeue the job in the case it fails.
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
conda activate trienv

srun python3 simulation_run_exp4.py -t 240215_experiment4_varying_factors_distance_boxheight