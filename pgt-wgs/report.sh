#!/bin/bash
set -euo pipefail

work_dir=$1
gene=$2
prefix=$3
pipeline=$(dirname $(readlink -f "$0"))

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
report=$pipeline/report.py

info=$work_dir/input.list
trans=$pipeline/etc/gene.bed
batch=$(basename $work_dir)

mkdir -p $work_dir/$batch
$python3 $report -info $info -trans $trans -work_dir $work_dir -gene $gene -prefix $prefix
cd $work_dir
tar -cf $batch.tar $batch
