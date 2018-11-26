#!/bin/bash
################################################################################
## partition: 
##   [12C] wmere, sbr, ivy, [16C] hwell, [8C] nehalem
##   Name    | Type | Nodes | Cores | Memory | Time limit | Others
##   --------|------|-------|-------|--------|------------|---------------------
##   k40     | GPU  | 24    | 20    | 64 GB  | 24 hour    | dual K40 GPUs
##   k20     | GPU  | 48    | 12    | 32 GB  | 24 hour    | single K20 GPU
##   hwell   | CPU  | 72    | 16    | 32 GB  | 48 hour    |
##   ivy     | CPU  | 80    | 12    | 16 GB  | 48 hour    |
##   sbr     | CPU  | 72    | 12    | 16 GB  | 48 hour    |
##   wmere   | CPU  | 96    | 12    | 12 GB  | 48 hour    |
##   nehalem | CPU  | 22    |  8    | 12 GB  | 48 hour    |
##
## qos (Quality of Service): 
##   veryshort ->  4 hour limit, increased priority
##   short     -> 12 hour limit, increased priority
##   long      ->  3  day limit, with a billing cost (UsageFactor)
##   medium    ->  1  day limit, increased priority
##   verylong  ->  7  day limit, with a billing cost (UsageFactor)
##   uberlong  -> 21  day limit, with a billing cost (UsageFactor)
################################################################################
#SBATCH -J TI
####SBATCH --time=20:00:00
#SBATCH --gres=gpu:p100
#SBATCH --partition=test
##SBATCH --constraints='x2680|x2695'
#SBATCH --nodes=1
##SBATCH --ntasks=112
##SBATCH --ntasks-per-core=1
##SBATCH --exclusive
#SBATCH --qos=veryshort # Max walltime for turbo is 8 hours
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

echo ""
echo "########################"
echo "#NN: Submit Dir.  ->" $SLURM_SUBMIT_DIR
echo "#NN: Job Nodelist ->" $SLURM_JOB_NODELIST
echo "#NN: No. of Nodes ->" $SLURM_JOB_NUM_NODES
echo "#NN: No. of Tasks ->" $SLURM_NTASKS
echo "########################"
echo ""

module purge
module load cuda/8.0.44
module load anaconda3


python umbrella.py

