#!/bin/bash
set -euo pipefail

pipeline=$1
work_dir=$2
gene=$3

export PATH=$pipeline/tools:$PATH
python3=$pipeline/tools/python3
vcftools=$pipeline/tools/vcftools
plink=$pipeline/tools/plink
tabix=$pipeline/tools/tabix
bgzip=$pipeline/tools/bgzip
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

trio=$(basename $work_dir)
anno_dir=$work_dir/Annotation
hap_dir=$work_dir/Hap
mkdir -p $anno_dir
mkdir -p $hap_dir

if [ -e $work_dir/shell/step13-Anno.sh.complete ]; then
    echo "Anno complete and skip"
else
    echo "`date` ADO Start"
    $python3 $ADO -vcf $work_dir/Family.SNP.INDEL.vcf.gz -target $target_gene -gene $gene -out $hap_dir/$trio.ADO.tsv -triosDepth 10 -embryoDepth 4
    echo "`date` ADO Done"

    echo "`date` Haplotype Start"
    $python3 $Haplotype -target $target_gene -gene $gene -vcf $work_dir/Family.SNP.INDEL.vcf.gz -correct neither -triosDepth 10 -embryoDepth 4 -prefix $hap_dir/$trio
    echo "`date` Haplotype Done"

    echo "`date` SnpSift Start"
    rm -f $hap_dir/gene.bed
    genes=(${gene//,/ })
    for g in ${genes[@]}
    do
        awk -v gene_name=$g '{if($1==gene_name) print $2"\t"$3"\t"$4}' $target_gene >> $hap_dir/gene.bed
    done
    $tabix -h -R $hap_dir/gene.bed $work_dir/Family.SNP.INDEL.vcf.gz > $hap_dir/Family.SNP.INDEL.filter.vcf
    $java -jar $snpEff/SnpSift.jar Annotate -tabix -noInfo $dbsnp $hap_dir/Family.SNP.INDEL.filter.vcf > $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf
    grep -v "##" $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf | cut -f1-5 > $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf.tsv
    $python3 $add_rsID $hap_dir/Family.SNP.INDEL.filter.dbsnp.vcf.tsv $hap_dir/$trio
    echo "`date` SnpSift Done"

    echo "`date` relationship Start"
    $vcftools --gzvcf $work_dir/Family.SNP.INDEL.vcf.gz --plink --out $work_dir/$trio
    $plink --noweb --file $work_dir/$trio --genome --out $work_dir/$trio
    $python3 $heatmap $work_dir/$trio.genome $trio $work_dir/$trio.relationship.pdf
    echo "`date` relationship Done"
    
    echo "`date` chromosome relationship Start"
    for chromosome in `cut -f1 $hap_dir/gene.bed | sort | uniq`
    do
        $tabix -h $work_dir/Family.SNP.INDEL.vcf.gz $chromosome | $bgzip > $work_dir/Family.SNP.INDEL.$chromosome.vcf.gz
        $vcftools --gzvcf $work_dir/Family.SNP.INDEL.$chromosome.vcf.gz --plink --out $work_dir/$trio.$chromosome
        $plink --noweb --file $work_dir/$trio.$chromosome --genome --out $work_dir/$trio.$chromosome
        $python3 $heatmap $work_dir/$trio.$chromosome.genome $trio.$chromosome $work_dir/$trio.$chromosome.relationship.pdf
    done
    echo "`date` chromosome relationship Done"

    echo "`date` snpEff Start"
    rm -f $anno_dir/gene.bed
    genes=(${gene//,/ })
    for g in ${genes[@]}
    do
        awk -v gene_name=$g '{if($6==gene_name) print $1"\t"$4"\t"$5}' $gene_trans | sort | uniq >> $anno_dir/gene.bed
    done
    $tabix -h -R $anno_dir/gene.bed $work_dir/Family.SNP.INDEL.vcf.gz > $anno_dir/Family.SNP.INDEL.filter.vcf
    $java -jar $snpEff/snpEff.jar -i vcf -o vcf -c $snpEff/snpEff.config hg19 $anno_dir/Family.SNP.INDEL.filter.vcf > $anno_dir/Family.SNP.INDEL.filter.ANN.vcf
    $python3 $Annotation -vcf $anno_dir/Family.SNP.INDEL.filter.ANN.vcf -out $anno_dir/Family.SNP.INDEL.filter.ANN.vcf.tsv
    $python3 $annotation_filter -trans $gene_trans -gene $gene -info $work_dir/../input.list -anno $anno_dir/Family.SNP.INDEL.filter.ANN.vcf.tsv -prefix $anno_dir/$trio -dipin $dipin
    echo "`date` snpEff Done"
fi
echo - | awk -v S=$SECONDS '{printf "step13-Anno\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' >> $work_dir/log
touch $work_dir/shell/step13-Anno.sh.complete
