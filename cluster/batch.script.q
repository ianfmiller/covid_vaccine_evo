#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 0-01:30:00
#SBATCH -J "covid.vacc"
#SBATCH --array=31-39
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu


INDEX=$(( $SLURM_ARRAY_TASK_ID + 0 ))
cd ~/Documents/GitHub/covid_vaccines_virulence_evolution/cluster/dir.$INDEX
srun R CMD BATCH analysis.cluster.R
