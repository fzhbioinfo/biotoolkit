#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
gatk=$pipeline/tools/gatk

mkdir -p $work_dir/javatmp
ls $work_dir/ChrVCF/*chr*.vcf.gz | awk '{print "-I "$0}' > $work_dir/vcf.list
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}.sh.complete

if [ -e $complete ]; then
    echo "gatk MergeVcfs complete and skip"
else
    echo "`date` gatk MergeVcfs Start"
    $gatk MergeVcfs \
    --arguments_file $work_dir/vcf.list -O $work_dir/Family.SNP.INDEL.vcf.gz --TMP_DIR $work_dir/javatmp
    echo "`date` gatk MergeVcfs Done"
fi
echo - | awk -v S=$SECONDS -v pipe=$pipe '{printf "%s\t%02d:%02d:%02d\n",pipe,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
