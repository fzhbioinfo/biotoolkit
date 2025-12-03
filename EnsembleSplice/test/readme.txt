python3 bgianno2vcf_splice.py ExprValid.uni.xlsx ExperimentValidateSplice.vcf
python3 ../maxentscan_wrapper.py -r /jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/NBS/hg19/hg19_chM_male_mask.fa -i ExperimentValidateSplice.vcf -p 2 --format_in vcf -o ExperimentValidateSplice.maxentscan.vcf
python3 ../scsnv_wrapper.py -i ExperimentValidateSplice.maxentscan.vcf --format_in vcf -p 2 -o ExperimentValidateSplice.maxentscan.scsnv.vcf
python3 ../spliceai_wrapper.py -i ExperimentValidateSplice.maxentscan.scsnv.vcf --format_in vcf -p 2 -o ExperimentValidateSplice.maxentscan.scsnv.spliceai.vcf
python3 parse_splice_result.py ExperimentValidateSplice.maxentscan.scsnv.spliceai.vcf ExperimentValidateSplice.maxentscan.scsnv.spliceai.vcf.result.tsv
