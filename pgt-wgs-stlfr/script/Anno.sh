#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
bedtools=$pipeline/tools/bedtools
vcftools=$pipeline/tools/vcftools
plink=$pipeline/tools/plink
tabix=$pipeline/tools/tabix
java=$pipeline/tools/java
Haplotype=$pipeline/stat/Haplotype.py
ADO=$pipeline/stat/ADO.py
heatmap=$pipeline/stat/heatmap.py
Annotation=$pipeline/stat/Annotation.py
annotation_filter=$pipeline/stat/annotation_filter.py
target_gene=$pipeline/etc/target.gene
gene_trans=$pipeline/etc/gene.bed
snpEff=$pipeline/snpEff
dbsnp=$pipeline/database/gatk/dbsnp_138.hg19.vcf.gz
add_rsID=$pipeline/stat/add_rsID.py
dipin=$pipeline/etc/dipin.common_name.txt

batch=$(basename $work_dir)
anno_dir=$work_dir/Annotation
hap_dir=$work_dir/Hap
mkdir -p $anno_dir
mkdir -p $hap_dir
info=$work_dir/input.list
genes=$work_dir/gene.list
script=$(basename $0)
pipe=`echo $script|sed 's/\.sh//g'`
complete=$work_dir/shell/${pipe}.sh.complete

if [ -e $complete ]; then
    echo "$pipe complete and skip"
else
    echo "`date` ADO Start"
    for g in `cat $genes`
    do
        cat $work_dir/LFR/*/Stat/$g.PS.bed | $bedtools sort -i - | $bedtools merge -i - > $hap_dir/$g.PS_merge.bed
    done
    $python3 $ADO -family_vcf $work_dir/Family.SNP.INDEL.vcf.gz -ps_bed_dir $hap_dir -gene $genes -out $hap_dir/$batch.ADO.tsv -triosDepth 10 -embryoDepth 4
    echo "`date` ADO Done"

    echo "`date` Haplotype Start"
    $python3 $Haplotype -info $info -lfr_dir $work_dir/LFR -target $target_gene -family_vcf $work_dir/Family.SNP.INDEL.vcf.gz -triosDepth 10 -embryoDepth 4 -prefix $hap_dir/$batch -extend
    echo "`date` Haplotype Done"

    echo "`date` SnpSift Start"
    cat $hap_dir/*.PS_merge.bed > $hap_dir/gene.PS.bed
    $tabix -h -R $hap_dir/gene.PS.bed $work_dir/Family.SNP.INDEL.vcf.gz > $hap_dir/Family.SNP.INDEL.filter.vcf
    $java -jar $snpEff/SnpSift.jar Annotate -tabix -noInfo $dbsnp $hap_dir/Family.SNP.INDEL.filter.vcf > $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf
    grep -v "##" $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf | cut -f1-5 > $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf.tsv
    $python3 $add_rsID $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf.tsv $hap_dir/$batch
    echo "`date` SnpSift Done"

    echo "`date` relationship Start"
    $vcftools --gzvcf $work_dir/Family.SNP.INDEL.vcf.gz --plink --out $work_dir/$batch
    $plink --noweb --file $work_dir/$batch --genome --out $work_dir/$batch
    $python3 $heatmap $work_dir/$batch.genome $batch $work_dir/$batch.relationship.pdf
    echo "`date` relationship Done"

    echo "`date` snpEff Start"
    rm -f $anno_dir/gene.bed
    for g in `cat $genes`
    do
        awk -v gene_name=$g '{if($6==gene_name) print $1"\t"$4"\t"$5}' $gene_trans | sort | uniq >> $anno_dir/gene.bed
    done
    $tabix -h -R $anno_dir/gene.bed $work_dir/Family.SNP.INDEL.vcf.gz > $anno_dir/Family.SNP.INDEL.filter.vcf
    $java -jar $snpEff/snpEff.jar -i vcf -o vcf -c $snpEff/snpEff.config hg19 $anno_dir/Family.SNP.INDEL.filter.vcf > $anno_dir/Family.SNP.INDEL.filter.ANN.vcf
    $python3 $Annotation -vcf $anno_dir/Family.SNP.INDEL.filter.ANN.vcf -out $anno_dir/Family.SNP.INDEL.filter.ANN.vcf.tsv
    $python3 $annotation_filter -trans $gene_trans -gene $genes -info $info -anno $anno_dir/Family.SNP.INDEL.filter.ANN.vcf.tsv -prefix $anno_dir/$batch -dipin $dipin
    echo "`date` snpEff Done"
fi
echo - | awk -v S=$SECONDS -v pipe=$pipe '{printf "%s\t%02d:%02d:%02d\n",pipe,S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $complete
