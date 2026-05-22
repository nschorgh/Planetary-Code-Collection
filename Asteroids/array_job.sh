#!/bin/bash

# Slurm script for array job with one integer as input argument

#SBATCH --job-name=myjob        # Job name
#SBATCH --output=stdout.%a      # Standard output and error log
#SBATCH --array=0-18
#SBATCH --ntasks=1              # Run a single task
#SBATCH --time=00:30:00         # Time limit hrs:min:sec

# --- Run Commands ---
a.out $SLURM_ARRAY_TASK_ID

## launch with 'sbatch array_job.sh'
