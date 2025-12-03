wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic_gaps.txt.gz
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
/jdfstj1/B2C_RD_P2/USR/fangzhonghai/software/anaconda3/envs/hgvs/bin/python genome_gap.py GCF_000001405.25_GRCh37.p13_genomic_gaps.txt.gz GRCh37 GRCh37.gap.regions.tsv
grep -v "#" GRCh37.gap.regions.tsv | cut -f1-3 | cat - GRCh37.Pseudoautosomal.regions.tsv | /jdfstj1/B2C_RD_P2/USR/fangzhonghai/software/local/bedtools sort -i - | /jdfstj1/B2C_RD_P2/USR/fangzhonghai/software/local/bedtools merge -d 1 -i - > GRCh37.N.regions.tsv
grep -v "@HD" ../database/hg19/ucsc.hg19.dict|cut -f2,3|sed -e 's/SN://g'|sed -e 's/LN:/1\t/g'|head -n 25 > GRCh37.regions.tsv
/jdfstj1/B2C_RD_P2/USR/fangzhonghai/software/local/bedtools subtract -a GRCh37.regions.tsv -b GRCh37.N.regions.tsv > GRCh37.bed
