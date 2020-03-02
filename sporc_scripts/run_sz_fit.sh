#!/bin/bash

jobname='SZ_Analysis'
jobfile='submit_catsrc.sh'
#cluster names need to be changed to accomodate multiple clusters at some point
name='rxj1347'
#this should be number of sims we want
nsims=100
#resolution and sgen need to be set before each run
resolution='nr'
sgen='3'

mkdir -p output

for nsim in $(seq 1 $nsims); do
    job=$jobname-$name-$sgen-$nsim-$resolution
    outfile=output/print_outputs_$name-$sgen-$nsim-$resolution.txt
    efile=output/error_outputs_$name-$sgen-$nsim-$resolution.txt
    export sgen;
    export name;
    export nsim;
    export resolution;
    sbatch --partition=tier3 -J $jobname -o $outfile -e $efile $jobfile
done

echo 'finished submitting jobs'