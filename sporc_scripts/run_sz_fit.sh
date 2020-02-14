#!/bin/bash 


jobname='SZ Cluster Analysis'
jobfile='submit_catsrc.sh'



#this is creating a list of cluster names
declare -a StringArray=('rjx1347' 'a0370')
#this is the number of sims that we want to iterate through
nsims=100
#resolution and sgen need to be set before each run
resolution='nr'
sgen='4'

mkdir -p output

for name in ${StringArray[@]}; do
    for nsim in $(seq 1 $lim_b); do
        job=$jobname-$name-$sgen$nsim
        outfile=output/output-a.$a-b.$b.txt
        export sgen
        export name
        export nsim
        export resolution
        sbatch --partition=debug -J $jobname -o $outfile $jobfile
    done;
done 

