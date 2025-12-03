#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
chromosome=$4

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk

chr_dir=$work_dir/$sample/chr
mkdir -p $work_dir/javatmp
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}-${chromosome}.sh.complete

if [ -e $complete ]; then
    echo "gatk MarkDuplicates $chromosome complete and skip"
else
    echo "`date` gatk MarkDuplicates $chromosome Start"
    $gatk MarkDuplicates \
    --CREATE_INDEX --CLEAR_DT false --TMP_DIR $work_dir/javatmp \
    -I $chr_dir/$chromosome.sort.bam -O $chr_dir/$chromosome.dup.bam -M $chr_dir/$chromosome.dup.metrics
    echo "`date` gatk MarkDuplicates $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v chromosome=$chromosome -v pipe=$pipe '{printf "%s-%s-%s\t%02d:%02d:%02d\n",pipe,sample,chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
