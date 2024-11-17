#!/bin/bash --login
#SBATCH --account=pawsey1103
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=200G
 
# Project directory
cd $MYSCRATCH/RDHB/RDHB_Occupancy

# Run setup
. src/cluster/a-setup-environment.sh
Rscript src/h-dynamic-spatial-occupancy_JAGS.R
