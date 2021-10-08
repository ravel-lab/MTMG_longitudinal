#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=32G -P jravel-lab -q all.q -N ZAPPS_virgo -j y -o /local/scratch/mfrance/logs/ -e /local/scratch/mfrance/logs/

#setting the number of jobs to be executed
#$ -t 1-1

cd /local/scratch/mfrance/PGATES/02_virgo/

infile=`sed -n "$SGE_TASK_ID p" MGs.lst`

rm temp_mapping/${infile}.out

/local/scratch/mfrance/PGATES/ZAPPS/02_virgo/runMapping.step1.sh -r /local/groupshare/ravel/mfrance/PGATE/ZAPPS/01_preprocess/2_fastq_se/${infile}.se.fq.gz -p $infile