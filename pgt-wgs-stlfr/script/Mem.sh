#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
fq1=$4
fq2=$5

export PATH=$pipeline/tools:$PATH
bwa=$pipeline/tools/bwa
samtools=$pipeline/tools/samtools
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa

filter_name=`echo $(basename $fq1) | sed 's/_1.fq.gz//g' `
filter_dir=$work_dir/$sample/Filter/$filter_name
aln_dir=$work_dir/$sample/Align
mkdir -p $aln_dir
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}-${filter_name}.sh.complete

if [ -e $complete ]; then
    echo "bwa mem $filter_name complete and skip"
else
    echo "`date` bwa mem $filter_name Start"
    $bwa mem \
    -M -t 12 -R "@RG\tID:PGD-MD\tSM:$sample\tLB:$sample\tPL:MGISEQ-2000" $genome \
    $filter_dir/${filter_name}_1.clean.fq.gz $filter_dir/${filter_name}_2.clean.fq.gz \
    | $samtools view -S -b -o $aln_dir/$filter_name.raw.bam -
    echo "`date` bwa mem $filter_name Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v filter_name=$filter_name -v pipe=$pipe '{printf "%s-%s-%s\t%02d:%02d:%02d\n",pipe,sample,filter_name,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
