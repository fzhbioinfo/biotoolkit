#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
merge_lfr_result=$pipeline/stat/merge_lfr_result.py

out_dir=$work_dir/LFR
mkdir -p $out_dir
info=$work_dir/input.list
batch=$(basename $work_dir)
mkdir -p $work_dir/$batch
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}.sh.complete

if [ -e $complete ]; then
    echo "$pipe complete and skip"
else
    echo "`date` merge_lfr_result Start"
    $python3 $merge_lfr_result -info $info -lfr_dir $out_dir -out $out_dir/$batch.stLFR.xlsx
    cp -f $out_dir/*/Stat/*.pdf $out_dir/$batch.stLFR.xlsx $work_dir/$batch
    cd $work_dir
    tar -cf $batch.tar $batch
    echo "`date` merge_lfr_result Done"
fi
echo - | awk -v S=$SECONDS -v pipe=$pipe '{printf "%s\t%02d:%02d:%02d\n",pipe,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
