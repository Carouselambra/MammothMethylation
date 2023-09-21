#!/bin/sh
#SBATCH -A snic2022-5-255
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -J Merge_Bin_elephants

cd /crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F
module load python3

python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_merge.py