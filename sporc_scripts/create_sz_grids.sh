#!/bin/bash

jobname='SZ_Spectral_Fitting'
jobfile='submit_grid.sh'

#input the name of the cluster
name='rxj1347'

#set the gridding for the lookup maps (will have dimension (N,N))
gridding=250
# create log file output directory above repo
mkdir -p /home/bjv7945/data/logs

#loop through the bands.
for value in BOLOCAM PMW PLW
do
    job=$jobname-$name-$value
    outfile=/home/bjv7945/data/logs/print_outputs_$name-$value.txt
    efile=/home/bjv7945/data/logs/error_outputs_$name-$value.txt
    export name;
    export value;
    export gridding;
    sbatch --partition=debug -J $jobname -o $outfile -e $efile $jobfile
done

echo 'finished submitting jobs'
