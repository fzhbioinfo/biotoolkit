#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
mut_type=$3

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa
hapmap=$pipeline/database/gatk/hapmap_3.3.hg19.sites.vcf.gz
omni=$pipeline/database/gatk/1000G_omni2.5.hg19.vcf.gz
snp=$pipeline/database/gatk/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
dbsnp=$pipeline/database/gatk/dbsnp_138.hg19.vcf.gz
indel=$pipeline/database/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

if [ $mut_type = 'SNP' ]; then
    if [ -e $work_dir/shell/step11-VariantRecalibrator-${mut_type}.sh.complete ]; then
        echo "gatk VariantRecalibrator complete and skip"
    else
        echo "`date` gatk VariantRecalibrator Start"
        $gatk SelectVariants \
        -select-type $mut_type --tmp-dir $work_dir/javatmp -R $genome \
        --variant $work_dir/temp.vcf.gz -O $work_dir/temp.${mut_type}.vcf
        $gatk VariantRecalibrator \
        -mode $mut_type --max-gaussians 4 --tmp-dir $work_dir/javatmp -R $genome \
        -V $work_dir/temp.${mut_type}.vcf -O $work_dir/temp.${mut_type}.recal.vcf \
        --tranches-file $work_dir/temp.${mut_type}.tranches \
        -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $snp \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp
        $gatk ApplyVQSR \
        -mode $mut_type --tmp-dir $work_dir/javatmp -R $genome --truth-sensitivity-filter-level 99.0 \
        -V $work_dir/temp.${mut_type}.vcf -O $work_dir/temp.${mut_type}.VQSR.vcf \
        --tranches-file $work_dir/temp.${mut_type}.tranches --recal-file $work_dir/temp.${mut_type}.recal.vcf
        $gatk SortVcf \
        --TMP_DIR $work_dir/javatmp -I $work_dir/temp.${mut_type}.VQSR.vcf -O $work_dir/temp.${mut_type}.VQSR.sort.vcf
        echo "`date` gatk VariantRecalibrator Done"
    fi
else
    if [ -e $work_dir/shell/step11-VariantRecalibrator-${mut_type}.sh.complete ]; then
        echo "gatk VariantRecalibrator complete and skip"
    else
        echo "`date` gatk VariantRecalibrator Start"
        $gatk SelectVariants \
        -select-type $mut_type --tmp-dir $work_dir/javatmp -R $genome \
        --variant $work_dir/temp.vcf.gz -O $work_dir/temp.${mut_type}.vcf
        $gatk VariantRecalibrator \
        -mode $mut_type --max-gaussians 4 --tmp-dir $work_dir/javatmp -R $genome \
        -V $work_dir/temp.${mut_type}.vcf -O $work_dir/temp.${mut_type}.recal.vcf \
        --tranches-file $work_dir/temp.${mut_type}.tranches \
        -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 $indel
        $gatk ApplyVQSR \
        -mode $mut_type --tmp-dir $work_dir/javatmp -R $genome --truth-sensitivity-filter-level 99.0 \
        -V $work_dir/temp.${mut_type}.vcf -O $work_dir/temp.${mut_type}.VQSR.vcf \
        --tranches-file $work_dir/temp.${mut_type}.tranches --recal-file $work_dir/temp.${mut_type}.recal.vcf
        $gatk SortVcf \
        --TMP_DIR $work_dir/javatmp -I $work_dir/temp.${mut_type}.VQSR.vcf -O $work_dir/temp.${mut_type}.VQSR.sort.vcf
        echo "`date` gatk VariantRecalibrator Done"
    fi
fi
echo - | awk -v S=$SECONDS -v mut_type=$mut_type '{printf "step11-VariantRecalibrator-%s\t%02d:%02d:%02d\n",mut_type,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step11-VariantRecalibrator-${mut_type}.sh.complete
