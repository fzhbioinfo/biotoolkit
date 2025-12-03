#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
chromosome=$3

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa
bed=$pipeline/etc/$chromosome.bed

vcf_dir=$work_dir/ChrVCF

if [ -e $work_dir/shell/step9-GenotypeGVCFs-${chromosome}.sh.complete ]; then
    echo "gatk GenotypeGVCFs $chromosome complete and skip"
else
    echo "`date` gatk GenotypeGVCFs $chromosome Start"
    $gatk GenotypeGVCFs \
    --tmp-dir $work_dir/javatmp -R $genome -L $bed \
    -V $vcf_dir/family.${chromosome}.gvcf.gz -O $vcf_dir/family.${chromosome}.vcf.gz
    echo "`date` gatk GenotypeGVCFs $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v chromosome=$chromosome '{printf "step9-GenotypeGVCFs-%s\t%02d:%02d:%02d\n",chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step9-GenotypeGVCFs-${chromosome}.sh.complete
