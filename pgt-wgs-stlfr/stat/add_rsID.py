# -*- coding:utf-8 -*-
import pandas as pd
import glob
import sys
import os
dbsnp = pd.read_csv(sys.argv[1], sep='\t')
dbsnp.rename(columns={'#CHROM': 'CHROM', 'ID': 'dbSNP'}, inplace=True)
columns = ['CHROM', 'POS', 'REF', 'ALT']
haps = glob.glob(f'{sys.argv[2]}*.hap.tsv')
for hap in haps:
    dir_name = os.path.dirname(hap)
    base_name = os.path.basename(hap).replace('hap.tsv', 'hap.xlsx')
    df = pd.read_csv(hap, sep='\t')
    df = pd.merge(df, dbsnp, on=columns, how='left')
    df.insert(2, 'dbSNP', df.pop('dbSNP'))
    df.to_excel(f'{dir_name}/{base_name}', index=False)
