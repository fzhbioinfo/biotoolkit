python3 /jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/EnsembleSplice/maxentscan_wrapper.py -r /jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/NBS/hg19/hg19_chM_male_mask.fa -i SPICE.c_to_g.0608.vcf -o SPICE.c_to_g.0608.MaxEntScan.vcf -p 2 --format_in vcf
python3 /jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/EnsembleSplice/maxentscan_wrapper.py -r /jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/NBS/hg19/hg19_chM_male_mask.fa -i SPICE.c_to_g.0608.vcf -o SPICE.c_to_g.0608.MaxEntScan.v2.vcf -p 2 --format_in vcf
python3 ../scsnv_wrapper.py -i SPICE.c_to_g.0608.MaxEntScan.v2.vcf -o SPICE.c_to_g.0608.MaxEntScan.v2.scSNV.vcf --format_in vcf -p 2
python3 ../spliceai_wrapper.py -i SPICE.c_to_g.0608.MaxEntScan.v2.scSNV.vcf -o SPICE.c_to_g.0608.MaxEntScan.v2.scSNV.SpliceAI.vcf --format_in vcf -p 2
python3 ../spliceai_wrapper.py -i SPICE.c_to_g.0608.MaxEntScan.v2.scSNV.vcf -o SPICE.c_to_g.0608.MaxEntScan.v2.scSNV.SpliceAI.vcf --format_in vcf -p 2
