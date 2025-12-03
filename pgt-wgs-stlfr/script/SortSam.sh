#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
fq1=$4
fq2=$5

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk

filter_name=`echo $(basename $fq1) | sed 's/_1.fq.gz//g' `
aln_dir=$work_dir/$sample/Align
mkdir -p $work_dir/javatmp
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}-${filter_name}.sh.complete

if [ -e $complete ]; then
    echo "gatk SortSam $filter_name complete and skip"
else
    echo "`date` gatk SortSam $filter_name Start"
    $gatk SortSam \
    -I $aln_dir/$filter_name.raw.bam -O $aln_dir/$filter_name.sort.bam \
    -SO coordinate --TMP_DIR $work_dir/javatmp --CREATE_INDEX
    echo "`date` gatk SortSam $filter_name Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v filter_name=$filter_name -v pipe=$pipe '{printf "%s-%s-%s\t%02d:%02d:%02d\n",pipe,sample,filter_name,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
