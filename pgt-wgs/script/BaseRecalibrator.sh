#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
chromosome=$4

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa
omni=$pipeline/database/gatk/1000G_omni2.5.hg19.vcf.gz
dbsnp=$pipeline/database/gatk/dbsnp_138.hg19.vcf.gz
snp=$pipeline/database/gatk/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
indel=$pipeline/database/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

chr_dir=$work_dir/$sample/chr

if [ -e $work_dir/shell/step6-BaseRecalibrator-${sample}-${chromosome}.sh.complete ]; then
    echo "gatk BaseRecalibrator $chromosome complete and skip"
else
    echo "`date` gatk BaseRecalibrator $chromosome Start"
    $gatk BaseRecalibrator \
    --tmp-dir $work_dir/javatmp -R $genome \
    -I $chr_dir/$chromosome.dup.bam -O $chr_dir/$chromosome.bqsr.recal_data.grp \
    --known-sites $omni --known-sites $dbsnp --known-sites $snp --known-sites $indel
    echo "`date` gatk BaseRecalibrator $chromosome Done"
    echo "`date` gatk ApplyBQSR $chromosome Start"
    $gatk ApplyBQSR \
    --tmp-dir $work_dir/javatmp -bqsr $chr_dir/$chromosome.bqsr.recal_data.grp \
    -I $chr_dir/$chromosome.dup.bam -O $chr_dir/$chromosome.bqsr.bam
    echo "`date` gatk ApplyBQSR $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v chromosome=$chromosome '{printf "step6-BaseRecalibrator-%s-%s\t%02d:%02d:%02d\n",sample,chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step6-BaseRecalibrator-${sample}-${chromosome}.sh.complete
