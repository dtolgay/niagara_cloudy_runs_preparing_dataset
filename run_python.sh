#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=cr_1_CO87_CII_H_O3_metallicity_above_minus_2
#SBATCH --output=cr_1_CO87_CII_H_O3_metallicity_above_minus_2.out
#SBATCH --error=cr_1_CO87_CII_H_O3_metallicity_above_minus_2.err

cd "/home/m/murray/dtolgay/scratch/cloudy_runs"

module purge 
ml python/3.11.5

python calculate_other_properties_from_finished_cloudy_runs.py