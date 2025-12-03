#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
sample=$3
chromosome=$4

export PATH=$pipeline/tools:$PATH
samtools=$pipeline/tools/samtools

aln_dir=$work_dir/$sample/Align
chr_dir=$work_dir/$sample/chr
mkdir -p $chr_dir
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}-${sample}-${chromosome}.sh.complete

if [ -e $complete ]; then
    echo "samtools merge $chromosome complete and skip"
else
    echo "`date` samtools merge $chromosome Start"
    ls $aln_dir/*.sort.bam > $chr_dir/$chromosome.sort.bam.list
    $samtools merge \
    -c -p -f -b $chr_dir/$chromosome.sort.bam.list -R $chromosome $chr_dir/$chromosome.sort.bam
    echo "`date` samtools merge $chromosome Done"
fi
echo - | awk -v S=$SECONDS -v sample=$sample -v chromosome=$chromosome -v pipe=$pipe '{printf "%s-%s-%s\t%02d:%02d:%02d\n",pipe,sample,chromosome,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
