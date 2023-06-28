#!/bin/bash
#SBATCH --partition=normal                      # "normal" partition
#SBATCH --nodes=1                               # we should only need one node for this!
#SBATCH --cpus-per-task=8                       # 16 CPUs will probably suffice.
#SBATCH --mem-per-cpu=16GB                      # we'll need an assload of memory for dealing with these huge matrices.
#SBATCH --output=./output/slurm/output/%x.out   # Output diagnostics and errors.
#SBATCH --export=./output/slurm/errors/%x.out
#SBATCH --time=5-00:00:00                       # Automatically set experiments to run for five days.
#SBATCH --mail-user=apizzime@gmu.edu            # Email me when things start/stop/error/get requeued.
#SBATCH --mail-type=ALL
#SBATCH --job-name=profile


module load gnu10
module load python
export PYTHONPATH=~/lib/python:$PYTHONPATH

sh profile.sh
