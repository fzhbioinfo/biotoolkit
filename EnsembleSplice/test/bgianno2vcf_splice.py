# -*- coding:utf-8 -*-
import pandas as pd
import pyfaidx
import sys


def bgi_anno_2_vcf_format(df, fa):
    df.reset_index(drop=True, inplace=True)
    df['#Chr'] = df['#Chr'].astype('str')
    if len(df[df['#Chr'].str.startswith('chr')]):
        df['#CHROM'] = df['#Chr']
    else:
        df['#CHROM'] = 'chr' + df['#Chr']
    df.loc[df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    df['ID'] = df['Splicing_effect']
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = '.'
    df['MuType'] = 'delins'
    df.loc[df['Ref'] == '.', 'MuType'] = 'ins'
    df.loc[df['Call'] == '.', 'MuType'] = 'del'
    df.loc[(df['Ref'].map(len) == 1) & (df['Call'].map(len) == 1) & (df['Ref'] != '.') & (df['Call'] != '.'), 'MuType'] = 'snp'
    df['POS'] = df['Stop']
    df.loc[df['MuType'] == 'del', 'POS'] = df.loc[df['MuType'] == 'del', 'Start']
    df.loc[df['MuType'] == 'delins', 'POS'] = df.loc[df['MuType'] == 'delins', 'Start']
    df['REF'] = df['Ref']
    df['ALT'] = df['Call']
    for i in range(df.shape[0]):
        if df.loc[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        elif df.loc[i, 'MuType'] == 'del':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'ALT'] = base
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
        elif df.loc[i, 'MuType'] == 'delins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        else:
            pass
    a = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].copy()
    a.sort_values(by=['#CHROM', 'POS'], ascending=True, inplace=True)
    b = df[['#Chr', 'Start', 'Stop', 'Ref', 'Call', '#CHROM', 'POS', 'REF', 'ALT']].copy()
    df.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO', 'MuType'], inplace=True)
    return a, b, df


def read_vcf_head(header):
    vh = open(header, 'r')
    fp_read = vh.read()
    vh.close()
    return fp_read


vcf_header = '/jdfstj1/B2C_RD_P2/USR/fangzhonghai/py_scripts/vcf_head.txt'
vcf_header_read = read_vcf_head(vcf_header)
hg19 = '/jdfstj1/B2C_RD_P2/USR/fangzhonghai/project/NBS/hg19/hg19_chM_male_mask.fa'
reference = pyfaidx.Fasta(hg19)

df_in = pd.read_excel(sys.argv[1])
#df_in.rename(columns={'Chr': '#Chr'}, inplace=True)
df_in = df_in[df_in['Ref'] != df_in['Call']].copy()

sample_vcf, _1, _2 = bgi_anno_2_vcf_format(df_in, reference)
with open(sys.argv[2], 'w') as f:
    f.write(vcf_header_read)
sample_vcf.to_csv(sys.argv[2], sep='\t', index=False, mode='a')
