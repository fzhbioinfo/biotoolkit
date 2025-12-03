#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
bcftools=$pipeline/tools/bcftools
tabix=$pipeline/tools/tabix

vcf_dir=$work_dir/ChrVCF

if [ -e $work_dir/shell/step10-Concat.sh.complete ]; then
    echo "bcftools concat complete and skip"
else
    echo "`date` bcftools concat Start"
    ls $vcf_dir/family.chr*.vcf.gz > $vcf_dir/vcf.list
    $bcftools concat \
    -a -D -q 30 -O z -f $vcf_dir/vcf.list -o $work_dir/temp.vcf.gz
    $tabix -f -p vcf $work_dir/temp.vcf.gz
    echo "`date` bcftools concat Done"
fi
echo - | awk -v S=$SECONDS '{printf "step10-Concat\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step10-Concat.sh.complete
