#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=96G -P jravel-lab -q threaded.q -pe thread 4 -N VCHIP_pre -j y -o /local/scratch/mfrance/logs/ -e /local/scratch/mfrance/logs/

#setting the number of jobs to be executed
#$ -t 1-1

cd /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess

infile=`sed -n "$SGE_TASK_ID p" MGs.lst`

#gunzip /local/scratch/mfrance/BIGEO/Irish/00_raw/human_removed/${infile}_R1.fastq.gz
#gunzip /local/scratch/mfrance/BIGEO/Irish/00_raw/human_removed/${infile}_R2.fastq.gz

./preprocessV2.sh $infile
#rm gzip /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/1_fastq_pe/${infile}*.gz
#rm gzip /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/2_fastq_se/${infile}*.gz

gzip /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/1_fastq_pe/${infile}*
gzip /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/2_fastq_se/${infile}*

mv /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/1_fastq_pe/${infile}* /local/groupshare/ravel/mfrance/PGATE/ZAPPS/01_preprocess/1_fastq_pe/
mv /local/scratch/mfrance/PGATES/ZAPPS/01_preprocess/2_fastq_se/${infile}* /local/groupshare/ravel/mfrance/PGATE/ZAPPS/01_preprocess/2_fastq_se/