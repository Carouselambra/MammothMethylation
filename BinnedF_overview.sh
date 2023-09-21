#!/bin/sh
#SBATCH -A snic2022-5-255
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 00:14:00
#SBATCH -J binf_overview_singles

cd /crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/binned_F
module load python3 

python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py

#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.merged.E.maximus.T13_T13_18 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.E.maximus.T13 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.E.maximus.T16
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.E.maximus.T18 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.L.africana.T11
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.americanum.ERR2260503 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.merged.M.primigenius.E467.L163.L164.P005.M6
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.primigenius.E467 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.primigenius.L163 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.primigenius.L164 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.primigenius.P005
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.M.primigenius.M6
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.merged.E.maximus.T13_T13_18 
#python3 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/BinnedF_overview.py RG.P.antiquus.ERR2260504
