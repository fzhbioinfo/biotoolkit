#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
chromosome=$4

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa
bed=$pipeline/etc/$chromosome.bed

chr_dir=$work_dir/$sample/chr
gvcf_dir=$work_dir/$sample/gvcf
mkdir -p $gvcf_dir

if [ -e $work_dir/shell/step7-HaplotypeCaller-${sample}-${chromosome}.sh.complete ]; then
    echo "gatk HaplotypeCaller $chromosome complete and skip"
else
    echo "`date` gatk HaplotypeCaller $chromosome Start"
    $gatk HaplotypeCaller \
    --tmp-dir $work_dir/javatmp -R $genome -ERC GVCF -L $bed \
    -I $chr_dir/$chromosome.bqsr.bam -O $gvcf_dir/$chromosome.gvcf.gz
    echo "`date` gatk HaplotypeCaller $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v chromosome=$chromosome '{printf "step7-HaplotypeCaller-%s-%s\t%02d:%02d:%02d\n",sample,chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step7-HaplotypeCaller-${sample}-${chromosome}.sh.complete
