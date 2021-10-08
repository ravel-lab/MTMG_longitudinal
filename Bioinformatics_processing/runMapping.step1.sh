#!/bin/bash

display_usage() {
	echo "Argument:"
	echo -e "\nUsage:\n-r read.fastq\n-p prefix \n"
	}
# if less than three arguments supplied, display usage
	if [  $# -le 1 ]
	then
		display_usage
		exit 1
	fi

# check whether user had supplied -h or --help . If yes display usage
	if [[ ( $# == "--help") ||  $# == "-h" ]]
	then
		display_usage
		exit 0
	fi

echo "Arguments are all good !!!"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--reads)
    READS="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--prefix)
    PREFIX="$2"
    shift # past argument
    shift # past value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}"

echo "reads file  = ${READS}"
echo "sample prefix  = ${PREFIX}"

dir=`echo $PWD`
output=$dir/temp_mapping
mkdir -p $output
reads=2
ref_db=/local/projects-t2/M8910/VOG_demo/0_db
dir_annot=/local/projects-t2/M8910/VOG_demo/1_VIRGO
ref=$ref_db/VIRGO
bowtie -p 16 -l 25 --fullref --chunkmbs 512 --best --strata -m 20 --mm $ref --suppress 2,4,5,6,7,8 $READS $output/$PREFIX.reads2ref
awk -F"\t" '{a[$2]++;}END{ for (i in a) print i"\t"a[i]}' $output/$PREFIX.reads2ref | sort -t$'\t' -k2nr,2 | awk -F"\t" -v reads=$reads '$2 >= reads' | awk -F"\t" 'NR==FNR{a[$1]=$2;next} ($1 in a) {print $0"\t"a[$1]}' $dir_annot/0.geneLength.txt - > $output/$PREFIX.out
rm $output/$PREFIX.reads2ref