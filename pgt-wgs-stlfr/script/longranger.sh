#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3

export PATH=$pipeline/tools:$PATH
longranger=$pipeline/longranger-2.2.2/longranger
reference=$pipeline/refdata-hg19-2.1.0
gatk=$pipeline/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar

out_dir=$work_dir/LFR
fq_dir=$out_dir/fq/$sample
mkdir -p $fq_dir
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}.sh.complete

if [ -e $complete ]; then
    echo "longranger $sample complete and skip"
else
    echo "`date` longranger $sample Start"
    cd $out_dir
    rm -rf $sample
    \time -v $longranger wgs --id=$sample --fastqs=$fq_dir --reference=$reference \
    --vcmode=gatk:${gatk} --localcores=96 --localmem=377
    echo "`date` longranger $sample Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v pipe=$pipe '{printf "%s-%s\t%02d:%02d:%02d\n",pipe,sample,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
