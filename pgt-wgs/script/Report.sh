#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
gene=$3

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
merge_QC=$pipeline/stat/merge_QC.py

hap_dir=$work_dir/Hap
anno_dir=$work_dir/Annotation
summary_dir=$work_dir/Summary
mkdir -p $summary_dir

if [ -e $work_dir/shell/step14-Report.sh.complete ]; then
    echo "Report complete and skip"
else
    $python3 $merge_QC $work_dir $work_dir/../input.list $gene
    echo "`date` Report Done"
    cp $work_dir/*.relationship.pdf $work_dir/*.QC.tsv $hap_dir/*.hap.xlsx $hap_dir/*.score.xlsx $anno_dir/*.anno.tsv $summary_dir
fi
echo - | awk -v S=$SECONDS '{printf "step14-Report\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step14-Report.sh.complete
