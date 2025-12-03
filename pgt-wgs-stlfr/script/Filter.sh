#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
fq1=$4
fq2=$5

export PATH=$pipeline/tools:$PATH
fastp=$pipeline/tools/fastp

filter_name=`echo $(basename $fq1) | sed 's/_1.fq.gz//g'`
filter_dir=$work_dir/$sample/Filter/$filter_name
mkdir -p $filter_dir
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}-${filter_name}.sh.complete

if [ -e $complete ]; then
    echo "fastp filter $filter_name complete and skip"
else
    echo "`date` fastp Filter $filter_name Start"
    $fastp \
    -i $fq1 -I $fq2 --thread 6 \
    -o $filter_dir/${filter_name}_1.clean.fq.gz -O $filter_dir/${filter_name}_2.clean.fq.gz \
    -j $filter_dir/$filter_name.fastp.json -h $filter_dir/$filter_name.fastp.html
    echo "`date` fastp filter $filter_name Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v filter_name=$filter_name -v pipe=$pipe '{printf "%s-%s-%s\t%02d:%02d:%02d\n",pipe,sample,filter_name,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
