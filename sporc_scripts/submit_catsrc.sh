#!/bin/bash -l

#SBATCH -A szeffect
#SBATCH --mail-user=vlb9398@g.rit.edu
##SBATCH --mail-type=ALL
#SBATCH -t 0-6:0:0
#SBATCH -p tier3 -n 1
#SBATCH --mem=600

cd ~/rsz
echo 'starting on:'
echo $name $sgen $nsim $resolution
python3 catsrc.py -run $name $sgen $nsim $resolution
