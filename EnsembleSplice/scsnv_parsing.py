# -*- coding:utf-8 -*-
import pandas as pd
import glob
import sys
import os
scSNV = glob.glob(sys.argv[1])
scSNV_df = pd.DataFrame()
work_dir = sys.argv[2]
chromosome = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
for snv in scSNV:
    df = pd.read_csv(snv, sep='\t', dtype={'chr': str, 'hg38_chr': str}, usecols=['chr', 'pos', 'ref', 'alt', 'hg38_chr', 'hg38_pos', 'ada_score', 'rf_score'], low_memory=False)
    df['chr'] = 'chr' + df['chr']
    df['hg38_chr'] = 'chr' + df['hg38_chr']
    scSNV_df = scSNV_df.append(df)
scSNV_df_hg19 = scSNV_df[['chr', 'pos', 'ref', 'alt', 'ada_score', 'rf_score']].copy()
scSNV_df_hg38 = scSNV_df[['hg38_chr', 'hg38_pos', 'ref', 'alt', 'ada_score', 'rf_score']].copy()
scSNV_df_hg19.rename(columns={'chr': '#chr'}, inplace=True)
scSNV_df_hg38.rename(columns={'hg38_chr': '#chr', 'hg38_pos': 'pos'}, inplace=True)
scSNV_df_hg19 = scSNV_df_hg19[scSNV_df_hg19['#chr'].isin(chromosome)].copy()
scSNV_df_hg38 = scSNV_df_hg38[scSNV_df_hg38['#chr'].isin(chromosome)].copy()
scSNV_df_hg19['pos'] = scSNV_df_hg19['pos'].astype('int64')
scSNV_df_hg38['pos'] = scSNV_df_hg38['pos'].astype('int64')
scSNV_df_hg19.sort_values(by=['#chr', 'pos'], ascending=[True, True], inplace=True)
scSNV_df_hg38.sort_values(by=['#chr', 'pos'], ascending=[True, True], inplace=True)
scSNV_df_hg19.to_csv(os.path.join(work_dir, 'dbscSNV1.1_hg19.tsv'), sep='\t', index=False)
scSNV_df_hg38.to_csv(os.path.join(work_dir, 'dbscSNV1.1_hg38.tsv'), sep='\t', index=False)
os.system('bgzip -f ' + os.path.join(work_dir, 'dbscSNV1.1_hg19.tsv'))
os.system('bgzip -f ' + os.path.join(work_dir, 'dbscSNV1.1_hg38.tsv'))
os.system('tabix -p vcf ' + os.path.join(work_dir, 'dbscSNV1.1_hg19.tsv.gz'))
os.system('tabix -p vcf ' + os.path.join(work_dir, 'dbscSNV1.1_hg38.tsv.gz'))
