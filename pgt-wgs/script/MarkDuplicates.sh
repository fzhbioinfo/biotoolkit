#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
chromosome=$4

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk

chr_dir=$work_dir/$sample/chr

if [ -e $work_dir/shell/step5-MarkDuplicates-${sample}-${chromosome}.sh.complete ]; then
    echo "gatk MarkDuplicates $chromosome complete and skip"
else
    echo "`date` gatk MarkDuplicates $chromosome Start"
    $gatk MarkDuplicates \
    --CREATE_INDEX --CLEAR_DT false --TMP_DIR $work_dir/javatmp \
    -I $chr_dir/$chromosome.sort.bam -O $chr_dir/$chromosome.dup.bam -M $chr_dir/$chromosome.dup.metrics
    echo "`date` gatk MarkDuplicates $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v chromosome=$chromosome '{printf "step5-MarkDuplicates-%s-%s\t%02d:%02d:%02d\n",sample,chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step5-MarkDuplicates-${sample}-${chromosome}.sh.complete
