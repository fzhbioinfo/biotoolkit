#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
bedtools=$pipeline/tools/bedtools
bed_dir=$pipeline/etc
Coverage=$pipeline/stat/Coverage.py
QC=$pipeline/stat/QC.py

chr_dir=$work_dir/$sample/chr
stat_dir=$work_dir/$sample/Stat
mkdir -p $stat_dir
genes=$work_dir/gene.list
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}.sh.complete

if [ -e $complete ]; then
    echo "$pipe complete and skip"
else
    echo "`date` Coverage wgs Start"
    $python3 $Coverage -bam_dir $chr_dir -bed_dir $bed_dir -process 5 -out $stat_dir/$sample.wgs.coverage.tsv
    echo "`date` Coverage wgs Done"

    echo "`date` Coverage gene Start"
    for g in `cat $genes`
    do
        cat $work_dir/LFR/*/Stat/$g.PS.bed | $bedtools sort -i - | $bedtools merge -i - > $stat_dir/$g.PS_merge.bed
        $python3 $Coverage -bam_dir $chr_dir -bed $stat_dir/$g.PS_merge.bed -out $stat_dir/$sample.$g.gene.coverage.tsv
    done
    echo "`date` Coverage gene Done"

    echo "`date` QC wgs Start"
    $python3 $QC -analysis_gender -coverage $stat_dir/$sample.wgs.coverage.tsv -out $stat_dir/$sample.wgs.qc.tsv
    echo "`date` QC wgs Done"

    echo "`date` QC gene Start"
    for g in `cat $genes`
    do
        $python3 $QC -coverage $stat_dir/$sample.$g.gene.coverage.tsv -out $stat_dir/$sample.$g.gene.qc.tsv
    done
    echo "`date` QC gene Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v pipe=$pipe '{printf "%s-%s\t%02d:%02d:%02d\n",pipe,sample,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
