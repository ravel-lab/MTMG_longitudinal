#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=16G -P jravel-lab -q all.q -N hr_gates -j y -o /local/scratch/mfrance/logs/ -e /local/scratch/mfrance/logs/

#setting the number of jobs to be executed
#$ -t 1-1

cd /local/scratch/mfrance/PGATES/ZAPPS/00_raw

infile=`sed -n -e "$SGE_TASK_ID p" MGs.lst`

raw_dir=/local/scratch/mfrance/PGATES/00_raw
human_genome_dir=/local/groupshare/ravel/mfrance/human_removal/source_hg38
out_dir=/local/projects-t3/PGATES/unzip_dump/

zcat $raw_dir/source/${infile}_R1.fastq.gz > $out_dir/${infile}_1.fastq
zcat $raw_dir/source/${infile}_R2.fastq.gz > $out_dir/${infile}_2.fastq

bmtagger.sh -b $human_genome_dir/GRch38_p12.bitmask -x $human_genome_dir/GRch38_p12.srprism -T /local/projects-t3/PGATES/hr/ -q 1 -1 $out_dir/${infile}_1.fastq -2 $out_dir/${infile}_2.fastq -X -o $raw_dir/human_removed/${infile}_de_human

echo $infile
wc -l $raw_dir/human_removed/${infile}_de_human_1.fastq