#!/bin/bash

#SBATCH --job-name=Array1d
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=08:00:00
#SBATCH --array=0-13

if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

things=$(cat sbatch_k2sc_all.txt)

if [ $SLURM_ARRAY_TASK_ID -lt ${#things[*]} ]; then
	module purge
	module load anaconda3/5.3.1

	virtualenv --system-site-packages py3.7.0
	source ~/pyenv/py3.7.0/bin/activate

    echo $SLURM_ARRAY_TASK_ID ${things[$SLURM_ARRAY_TASK_ID]}
fi