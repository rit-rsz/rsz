#!/bin/bash -l

#SBATCH -A szeffect
#SBATCH --mail-user=vlb9398@g.rit.edu
##SBATCH --mail-type=ALL
#SBATCH -t 12:0:0
#SBATCH -p tier3 -n 1
#SBATCH --mem=500M

cd ~/rsz/sz/
python3 clus_sz_grids.py -run $value $name $gridding
