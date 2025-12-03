#!/bin/bash
set -euo pipefail

pro_dir=$(dirname $0)
highconf_bed=$pro_dir/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.add_chr.bed
benchmark_vcf=$pro_dir/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.add_chr.vcf.gz
standard_stat=$pro_dir/standard_stat.py

target_bed=$1
vcf_in=$2
prefix=$3
sample=$(basename $prefix)

bedtools intersect -a $highconf_bed -b $target_bed | sort -k 1,1V -k 2,2n -k 3,3n | bedtools merge -i - > $prefix.highconf_target.bed
tabix -h -R $prefix.highconf_target.bed $benchmark_vcf | bgzip -f > $prefix.highconf_target.benchmark.vcf.gz
tabix -f -p vcf $prefix.highconf_target.benchmark.vcf.gz
tabix -h -R $prefix.highconf_target.bed $vcf_in | bgzip -f > $prefix.highconf_target.input.vcf.gz
tabix -f -p vcf $prefix.highconf_target.input.vcf.gz
gatk Concordance -eval $prefix.highconf_target.input.vcf.gz --truth $prefix.highconf_target.benchmark.vcf.gz --summary $prefix.Concordance.txt -tpfn $prefix.tpfn.vcf.gz -tpfp $prefix.tpfp.vcf.gz
python3 $standard_stat $prefix.Concordance.txt $sample $prefix.highconf_target.bed $prefix.standard_stat.json
