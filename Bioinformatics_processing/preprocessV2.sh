#!/bin/bash
## qsub -cwd -b y -l mem_free=32G -P jravel-lab -q all.q -N pre -V -j y -o logs/ ./preprocess_MGMT.sh
## qsub -cwd -b y -l mem_free=4G -P jravel-lab -q threaded.q@galactus -V -j y -o logs/ ./B_QCFirst

########do human removal first for vaginal sample ######
#mkdir -p 0_humanRemove
############RUN BMTAGGER !!!!

sample=$1
miniLen=50 ### 151x0.5

#human removal result dir
dir=/local/scratch/mfrance/PGATES/00_raw/human_removed/
tmpdir=_temp; mkdir -p $tmpdir
sedir=2_fastq_se; mkdir -p $sedir
pedir=1_fastq_pe; mkdir -p $pedir
#fadir=3_fa; mkdir -p $fadir
#mkdir -p 02_fasta/
	infile1=$dir/${sample}_de_human_1.fastq
	#tmpfile1=${sample}_1_tmp.fastq
	infile2=$dir/${sample}_de_human_2.fastq
        #gunzip $infile1
        #gunzip $infile2
	#tmpfile2=${sample}_2_tmp.fastq
	#### remove ribosomal RNA ######
	/local/projects-t2/M8910/software/sortmerna-2.1b/scripts/merge-paired-reads.sh $infile1 $infile2 $tmpdir/$sample.interleaved.fq
	/home/bma/software2/sortmerna-2.1b/sortmerna --ref /home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-bac-16s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-bac-23s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-arc-16s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-arc-23s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-euk-18s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/home/bma/software2/sortmerna-2.1b/index/silva-euk-28s:/home/bma/software2/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/home/bma/software2/sortmerna-2.1b/index/rfam-5s-db:/home/bma/software2/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/home/bma/software2/sortmerna-2.1b/index/rfam-5.8s-db --reads $tmpdir/$sample.interleaved.fq --num_alignments 1 --fastx --aligned $tmpdir/$sample.reads_rRNA --other $tmpdir/$sample.reads_non_rRNA --log -a 8 -m 2400 --paired_in -v --sam
	/local/projects-t2/M8910/software/sortmerna-2.1b/scripts/unmerge-paired-reads.sh $tmpdir/$sample.reads_non_rRNA.fq $tmpdir/$sample.reads_non_rRNA.R1.fq $tmpdir/$sample.reads_non_rRNA.R2.fq
	#cat $tmpdir/$sample.reads_non_rRNA.fq | /usr/local/bin/fastq-stats > $tmpdir/$sample.nonRibo.stat #count one given paired reads
	rm $tmpdir/$sample.reads_rRNA.sam $tmpdir/$sample.reads_non_rRNA.fq $tmpdir/$sample.interleaved.fq $tmpdir/$sample.reads_rRNA.fq
	#cat $X | grep "^reads" | awk -F"\t" '{print $2}' | awk -v name=$name '{sum=+$1}END{print name"\t"sum}'; done > PE1.stat.txt

	#### remove phiX reads ####
        #bowtie -p 4 -l 21 --mm /local/groupshare/ravel/BMA/1_db/phiX/phiX --suppress 1,2,3,4,5,6,7,8 $tmpdir/${sample}.reads_non_rRNA.R1.fq --un $tmpdir/${sample}.reads_non_rRNA.R1.fq
        #bowtie -p 4 -l 21 --mm /local/groupshare/ravel/BMA/1_db/phiX/phiX --suppress 1,2,3,4,5,6,7,8 $tmpdir/${sample}.reads_non_rRNA.R2.fq --un $tmpdir/${sample}.reads_non_rRNA.R2.fq

        #cat $tmpdir/${sample}.reads_non_rRNA.R1.fq | /usr/local/bin/fastq-stats/fastq-stats > $tmpdir/$sample.phiX.stat
	### quality control #######
	###turns out trimmonic PE mode does not generate paired end reads in the same order, not even the same set of reads!!!!
	java -jar /home/mfrance/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 $tmpdir/$sample.reads_non_rRNA.R1.fq $tmpdir/$sample.reads_non_rRNA.R2.fq $pedir/${sample}_1_de_all.fastq $pedir/${sample}_unpaired.1.fastq $pedir/${sample}_2_de_all.fastq $pedir/${sample}_unpaired.2.fastq ILLUMINACLIP:/home/mfrance/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:6:10 MINLEN:$miniLen
	rm $tmpdir/$sample.reads_non_rRNA.R1.fq $tmpdir/$sample.reads_non_rRNA.R2.fq

	#fastp
    #/home/mfrance/software/fastp/fastp -w 4 --html $tmpdir/${sample}_report.html --json $tmpdir/${sample}_report.json -g -l $miniLen -W 4 -M 20 -i $tmpdir/$sample.reads_non_rRNA.R1.fq -I $tmpdir/$sample.reads_non_rRNA.R2.fq -o $pedir/${sample}_1_de_all.fastq -O $pedir/${sample}_2_de_all.fastq --unpaired1 $pedir/${sample}_unpaired.1.fastq --unpaired2 $pedir/${sample}_unpaired.2.fastq
    #rm $tmpdir/$sample.reads_non_rRNA.R1.fq $tmpdir/$sample.reads_non_rRNA.R2.fq
        ###nead to merge the trimmed fq files together to form one file
	cat $pedir/${sample}_1_de_all.fastq $pedir/${sample}_unpaired.1.fastq $pedir/${sample}_2_de_all.fastq $pedir/${sample}_unpaired.2.fastq > $sedir/${sample}.se.fq
	cat $pedir/${sample}_unpaired.1.fastq $pedir/${sample}_unpaired.2.fastq > $pedir/${sample}.unpaired.fq
	cat $sedir/${sample}.se.fq | /usr/local/bin/fastq-stats > $sedir/${sample}.total.stat

	###### to match up paired end reads and sort in the same order #######
	/home/mfrance/software/MiSeq16S/fq_getPairAndOrphan1.8.py $pedir/${sample}_1_de_all.fastq $pedir/${sample}_2_de_all.fastq $pedir/${sample}_de_all_1.fastq $pedir/${sample}_de_all_2.fastq $pedir/${sample}_orphan.fastq

	####fq_getPairAndOrphan1.8 does not generate the same order always, while it does generate same set of reads
	cat $pedir/${sample}_de_all_1.fastq | /home/mfrance/software/MiSeq16S/fq_mergelines.pl > $pedir/$sample.R1.tab
	cat $pedir/${sample}_de_all_2.fastq | /home/mfrance/software/MiSeq16S/fq_mergelines.pl > $pedir/$sample.R2.tab
	sort -k1,1 $pedir/$sample.R1.tab -T /local/groupshare/ravel/mfrance/temp | /home/mfrance/software/MiSeq16S/fq_splitlines.pl > $pedir/${sample}.R1.fq
	sort -k1,1 $pedir/$sample.R2.tab -T /local/groupshare/ravel/mfrance/temp | /home/mfrance/software/MiSeq16S/fq_splitlines.pl > $pedir/${sample}.R2.fq
	rm $pedir/$sample.R1.tab $pedir/$sample.R2.tab

        ##### assemble PE for assembly ####
	#fastq-join -m 30 $pedir/${sample}.R1.fq $pedir/${sample}.R2.fq -o $sedir/${sample}.%.fq
	#fq2fa $sedir/${sample}.join.fq $fadir/$sample.join.fa
	rm $sedir/$sample.un1.fq $sedir/$sample.un2.fq

	###### generate one single fasta file ######
	#/usr/local/bin/fq2fa $sedir/${sample}_se.fq $fadir/$sample.fa
	#interleaved=$fadir/$sample.interleaved.fa
        #fq2fa --filter --merge $pedir/${sample}.R1.fq $pedir/${sample}.R2.fq $interleaved
	#seqtk seq -a $qcdir/${sample}_de_all.shuffle.fq > 02_fasta/${sample}_de_all.fa

rm $pedir/${sample}_1_de_all.fastq $pedir/${sample}_2_de_all.fastq $pedir/${sample}_orphan.fastq $pedir/${sample}_unpaired.?.fastq $pedir/${sample}_de_all_?.fastq #$tmpdir/$tmpfile1 $tmpdir/$tmpfile2 $tmpdir/${sample}_1_tmp?.fastq
