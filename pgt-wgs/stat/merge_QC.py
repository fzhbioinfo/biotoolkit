# -*- coding:utf-8 -*-
from collections import defaultdict
import pandas as pd
import glob
import sys
import os
work_dir = sys.argv[1]
info = pd.read_csv(sys.argv[2], sep='\t', usecols=['sample', 'sampleID', 'relation', 'gender'])
genes = sys.argv[3]
info.rename(columns={'gender': 'gender_provided'}, inplace=True)
info.drop_duplicates(subset=['sampleID'], inplace=True)
info['ID'] = info['sample'] + '_' + info['sampleID']
trio = os.path.basename(work_dir)
info = info.loc[(info['sample'].isin(['Father', 'Mother'])) | (info['ID'] == trio) | (info['sample'].str.contains('Embryo'))]
# 合并wgs范围和基因qc
qc_wgs = glob.glob(f'{work_dir}/*/Stat/*.wgs.qc.tsv')
qc_wgs_df = pd.DataFrame()
qc_gene_list = list()
for qc in qc_wgs:
    sample = os.path.basename(qc).split('.')[0]
    df = pd.read_csv(qc, sep='\t')
    df['sample'] = sample
    qc_wgs_df = qc_wgs_df.append(df, sort=False)
    qc_gene_dic = {'sample': sample}
    for gene in genes.split(','):
        qc_gene = pd.read_csv(f'{work_dir}/{sample}/Stat/{sample}.{gene}.gene.qc.tsv', sep='\t')
        qc_gene_dic[f'mean_depth(rm_dup)_{gene}'] = qc_gene.loc[0, 'mean_depth(rm_dup)']
    qc_gene_list.append(qc_gene_dic)
qc_wgs_df = pd.merge(info, qc_wgs_df, on=['sample'], how='left')
qc_gene_df = pd.DataFrame(qc_gene_list)
qc_df = pd.merge(qc_wgs_df, qc_gene_df, on=['sample'], how='left')
qc_df.fillna('.', inplace=True)
# 合并ADO
ado_dic = defaultdict(dict)
ado = pd.read_csv(f'{work_dir}/Hap/{trio}.ADO.tsv', sep='\t', index_col=0)
for row in ado.itertuples():
    sample = row.Index
    if sample in ['Father', 'Mother']:
        continue
    ado_dic[sample]['sample'] = sample
    try:
        rate = round(int(row.ADO) / int(row.meet), 4)
    except ZeroDivisionError:
        rate = '.'
    ado_dic[sample][f'ADO_{row.gene}'] = rate
ado_df = pd.DataFrame.from_dict(ado_dic, orient='index')
if ado_df.empty:
    qc_df['ADO_all_gene'] = '.'
    for gene in genes.split(','):
        qc_df[f'ADO_{gene}'] = '.'
    out_df = qc_df
else:
    out_df = pd.merge(qc_df, ado_df, on=['sample'], how='left')
    out_df.fillna('.', inplace=True)
out_df.to_csv(f'{work_dir}/{trio}.QC.tsv', sep='\t', index=False)
