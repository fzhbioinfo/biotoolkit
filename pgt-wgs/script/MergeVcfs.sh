#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk

snp=$work_dir/temp.SNP.VQSR.sort.vcf
indel=$work_dir/temp.INDEL.VQSR.sort.vcf

if [ -e $work_dir/shell/step12-MergeVcfs.sh.complete ]; then
    echo "gatk MergeVcfs complete and skip"
else
    echo "`date` gatk MergeVcfs Start"
    $gatk MergeVcfs \
    -I $snp -I $indel -O $work_dir/Family.SNP.INDEL.vcf.gz --TMP_DIR $work_dir/javatmp
    echo "`date` gatk MergeVcfs Done"
fi
echo - | awk -v S=$SECONDS '{printf "step12-MergeVcfs\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step12-MergeVcfs.sh.complete
