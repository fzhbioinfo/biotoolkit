#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
gene=$4

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
bed_dir=$pipeline/etc
target=$pipeline/etc/target.gene
Coverage=$pipeline/stat/Coverage.py
QC=$pipeline/stat/QC.py
gene_stat=$pipeline/stat/gene_stat.py
genome=$pipeline/database/hg19/hg19_chM_male_mask.fa

chr_dir=$work_dir/$sample/chr
stat_dir=$work_dir/$sample/Stat
mkdir -p $stat_dir
genes=(${gene//,/ })

if [ -e $work_dir/shell/step7-Stat-${sample}.sh.complete ]; then
    echo "Coverage stat complete and skip"
else
    echo "`date` Coverage wgs Start"
    $python3 $Coverage -bam_dir $chr_dir -bed_dir $bed_dir -process 5 -out $stat_dir/$sample.wgs.coverage.tsv
    echo "`date` Coverage wgs Done"
    
    echo "`date` Coverage gene Start"
    for g in ${genes[@]}
    do
        awk -v gene_name=$g '{if($1==gene_name) print $2"\t"$3"\t"$4}' $target > $stat_dir/$g.bed
        $python3 $Coverage -bam_dir $chr_dir -bed $stat_dir/$g.bed -out $stat_dir/$sample.$g.gene.coverage.tsv
    done
    echo "`date` Coverage gene Done"

    echo "`date` QC wgs Start"
    $python3 $QC -analysis_gender -coverage $stat_dir/$sample.wgs.coverage.tsv -out $stat_dir/$sample.wgs.qc.tsv
    echo "`date` QC wgs Done"
    
    echo "`date` QC gene Start"
    for g in ${genes[@]}
    do
        $python3 $gene_stat $stat_dir/$sample.$g.gene.coverage.tsv $genome $stat_dir/$sample.$g.gene.qc.tsv
    done
    echo "`date` QC gene Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample '{printf "step7-Stat-%s\t%02d:%02d:%02d\n",sample,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step7-Stat-${sample}.sh.complete
