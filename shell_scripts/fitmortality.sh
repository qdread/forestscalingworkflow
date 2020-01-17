#!/bin/bash --login
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=NONE
#SBATCH --array=1-3

# Hardcoded all options.
# chain is the slurm array task id

((seed=777+SLURM_ARRAY_TASK_ID))

~/forestlight/stancode/mortreg_fg_v3 \
  sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 \
  random seed=${seed} data file=~/forestlight/stanrdump/mortalitydump.r \
  output file=~/forestlight/stanoutput/fit_mortality_${SLURM_ARRAY_TASK_ID}.csv

