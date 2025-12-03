#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
bed_dir=$pipeline/etc
target=$pipeline/etc/target.gene
region=$pipeline/etc/region.tsv
Coverage=$pipeline/stat/Coverage.py
ps_coordinate=$pipeline/stat/ps_coordinate.py
QC=$pipeline/stat/QC.py
ps_stat=$pipeline/stat/ps_stat.py
ps_plot=$pipeline/stat/ps_plot.py
ps_detect=$pipeline/stat/ps_detect.py

out_dir=$work_dir/LFR
stat_dir=$out_dir/$sample/Stat
mkdir -p $stat_dir
bam=$out_dir/$sample/outs/phased_possorted_bam.bam
vcf=$out_dir/$sample/outs/phased_variants.vcf.gz
info=$work_dir/input.list
genes=$work_dir/gene.list
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}.sh.complete

if [ -e $complete ]; then
    echo "$pipe complete and skip"
else
    echo "`date` Coverage wgs Start"
    $python3 $Coverage -bam $bam -bed_dir $bed_dir -process 5 -out $stat_dir/$sample.wgs.coverage.tsv
    echo "`date` Coverage wgs Done"

    echo "`date` Coverage gene Start"
    for g in `cat $genes`
    do
        $python3 $ps_coordinate -target $target -gene $g -lfr $vcf -out $stat_dir/$g.PS.bed
        $python3 $Coverage -bam $bam -bed $stat_dir/$g.PS.bed -out $stat_dir/$sample.$g.gene.coverage.tsv
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

    echo "`date` ps_stat Start"
    $python3 $ps_stat -lfr $vcf -region $region -out $stat_dir/$sample.PS_stat.tsv
    echo "`date` ps_stat Done"

    echo "`date` ps_plot Start"
    $python3 $ps_plot -lfr $vcf -target $target -gene $genes -out $stat_dir/$sample.PS_plot.pdf
    echo "`date` ps_plot Done"

    echo "`date` ps_detect Start"
    $python3 $ps_detect -info $info -sample $sample -lfr_workdir $out_dir/$sample -out $stat_dir/$sample.PS_detect.tsv
    echo "`date` ps_detect Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v pipe=$pipe '{printf "%s-%s\t%02d:%02d:%02d\n",pipe,sample,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
