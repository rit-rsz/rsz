#!/bin/bash

jobname='Lense_tmplt'
jobfile='submit_lense.sh'
name='rxj1347'
chunks=$1
band=$2

mkdir -p output

for chunk in $(seq 1 $chunks); do
    job=$jobname-$name-$chunk-$band
    outfile=output/print_$name-$band-$chunk.txt
    efile=output/error_$name-$band-$chunk.txt
    export name;
    export chunk;
    export band;
    sbatch --partition=tier3 -J $jobname -o $outfile -e $efile $jobfile
done

echo 'finished submitting jobs'
