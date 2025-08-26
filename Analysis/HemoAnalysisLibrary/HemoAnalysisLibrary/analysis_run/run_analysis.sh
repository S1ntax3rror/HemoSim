#!/bin/bash
#SBATCH --job-name=Analysis_Main
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=2000

#gpu05
echo $HOSTNAME

python3 analysis_runner.py

