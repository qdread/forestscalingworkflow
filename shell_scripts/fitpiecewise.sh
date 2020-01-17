#!/bin/bash --login
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=NONE
#SBATCH --array=1-3

# Include model as variable (density1, production1, etc.) with --export
# Include scaling, guild and year as variables with --export
# Include sampling (NS) and warmup (NW) numbers with --export
# Include (random) seed with --export
# chain is the slurm array task id

((seed=seed+SLURM_ARRAY_TASK_ID))

~/forestlight/stancode/${model} \
  sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 random seed=${seed} \
  data file=~/forestlight/stanrdump/dump_${scaling}_${guild}_${year}.r \
  output file=~/forestlight/stanoutput/fit_${model}_${scaling}_${guild}_${year}_${SLURM_ARRAY_TASK_ID}.csv
