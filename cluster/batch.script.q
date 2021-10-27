#!/bin/bash
#
#SBATCH -N 1
#SBATCH --cores-per-socket=20
#SBATCH -J "covid.vacc"
#SBATCH --array=1-12
#SBATCH --mem=20g
#SBATCH -t 0-05:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu


INDEX=$(( $SLURM_ARRAY_TASK_ID + 0 ))
cd ~/Documents/GitHub/covid_vaccines_virulence_evolution/cluster/dir.$INDEX
srun R CMD BATCH analysis.cluster.R
