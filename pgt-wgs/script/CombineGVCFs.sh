#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
chromosome=$3

export PATH=$pipeline/tools:$PATH
java=$pipeline/tools/java
GenomeAnalysisTK=$pipeline/tools/GenomeAnalysisTK.jar
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa
bed=$pipeline/etc/$chromosome.bed

vcf_dir=$work_dir/ChrVCF
mkdir -p $vcf_dir

if [ -e $work_dir/shell/step8-CombineGVCFs-${chromosome}.sh.complete ]; then
    echo "gatk CombineGVCFs $chromosome complete and skip"
else
    echo "`date` gatk CombineGVCFs $chromosome Start"
    ls $work_dir/*/gvcf/$chromosome.gvcf.gz | awk '{print "-V "$0}' > $vcf_dir/$chromosome.gvcf.list
    $java -jar $GenomeAnalysisTK -T CombineGVCFs \
    -L $bed -R $genome -l ERROR \
    -args $vcf_dir/$chromosome.gvcf.list -o $vcf_dir/family.${chromosome}.gvcf.gz
    echo "`date` gatk CombineGVCFs $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v chromosome=$chromosome '{printf "step8-CombineGVCFs-%s\t%02d:%02d:%02d\n",chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step8-CombineGVCFs-${chromosome}.sh.complete
