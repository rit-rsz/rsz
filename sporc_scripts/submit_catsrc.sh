#!/bin/bash -l

#SBATCH -A szeffect
#SBATCH --mail-user bjv7945@rit.edu
#SBATCH --mail-type=ALL
#SBATCH -t 0-3:0:5
#SBATCH -p debug -n 1	
#SBATCH --mem=600
#SBATCH -o output/output-$name-$sgen$nsim.txt
#SBATCH -e output/error-$name-$sgen$nsim.txt

cd rsz/
python3 catsrc -run $name $sgen $nsim $resolution 
