#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
fq1=$4
fq2=$5

export PATH=$pipeline/tools:$PATH
splitBarcode=$pipeline/splitBarcode/split10x
barcode=$pipeline/splitBarcode/etc/barcode.list
map=$pipeline/splitBarcode/etc/4M-with-alts-february-2016.txt.gz

out_dir=$work_dir/LFR
fq_dir=$out_dir/fq/$sample
mkdir -p $fq_dir
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}.sh.complete

if [ -e $complete ]; then
    echo "splitBarcode $sample complete and skip"
else
    echo "`date` splitBarcode $sample Start"
    \time -v $splitBarcode -fq1 $fq1 -fq2 $fq2 -prefix $fq_dir/$sample \
    -bc $barcode -map $map
    echo "`date` splitBarcode $sample Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v pipe=$pipe '{printf "%s-%s\t%02d:%02d:%02d\n",pipe,sample,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
