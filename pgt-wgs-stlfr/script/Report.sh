#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
merge_QC=$pipeline/stat/merge_QC.py
report=$pipeline/stat/report.py

info=$work_dir/input.list
batch=$(basename $work_dir)
mkdir -p $work_dir/$batch
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}.sh.complete

if [ -e $complete ]; then
    echo "$pipe complete and skip"
else
    echo "`date` merge_QC Start"
    $python3 $merge_QC -work_dir $work_dir -info $info -out $work_dir/$batch.QC.xlsx
    echo "`date` merge_QC Done"

    echo "`date` report Start"
    $python3 $report -work_dir $work_dir -info $info -out $work_dir/$batch.report.xlsx -lfr_qc $work_dir/LFR/$batch.stLFR.xlsx -wgs_qc $work_dir/$batch.QC.xlsx 
    cp $work_dir/*.xlsx $work_dir/*.pdf $work_dir/Hap/*.hap.* $work_dir/$batch
    rm -f $work_dir/$batch/*.tsv
    cd $work_dir
    tar -cf $batch.tar $batch
    echo "`date` report Done"
fi
echo - | awk -v S=$SECONDS -v pipe=$pipe '{printf "%s\t%02d:%02d:%02d\n",pipe,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
