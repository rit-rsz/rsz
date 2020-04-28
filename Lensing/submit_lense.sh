#!/bin/bash -l

#SBATCH -A szeffect
#SBATCH --mail-user=vlb9398@g.rit.edu
##SBATCH --mail-type=ALL
#SBATCH -t 0-0:30:0
#SBATCH -p tier3 -n 1
#SBATCH --mem=1G

cd ~/rsz/Lensing
echo 'starting on:'
echo $name $chunk $band
python3 lense_temp_step2.py -run $name $chunk $band
