#!/bin/bash -l

#SBATCH -A szeffect
#SBATCH --mail-user bjv7945@rit.edu
#SBATCH --mail-type=ALL
#SBATCH -t 0-0:0:5
#SBATCH -p debug -n 1	
#SBATCH --mem=600


cd ~/rsz
echo 'starting on:' 
echo $name $sgen-$nsim $resolution
python3 catsrc.py $name $sgen $nsim $resolution
