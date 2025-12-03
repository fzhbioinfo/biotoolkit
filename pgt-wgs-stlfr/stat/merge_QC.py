# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
import glob
import os


def merge_qc(work_dir, samples):
    # merge wgs qc
    df_wgs_qc = pd.DataFrame()
    for wgs_qc in glob.glob(f'{work_dir}/*/Stat/*.wgs.qc.tsv'):
        tmp = pd.read_csv(wgs_qc, sep='\t')
        sample = os.path.basename(os.path.dirname(os.path.dirname(wgs_qc)))
        tmp['sample'] = sample
        df_wgs_qc = df_wgs_qc.append(tmp)
    # merge gene qc
    df_gene_qc = pd.DataFrame()
    for sample in samples:
        tmp_ = pd.DataFrame()
        for gene_qc in glob.glob(f'{work_dir}/{sample}/Stat/*.gene.qc.tsv'):
            cols = ['mean_depth(rm_dup)', 'coverage>=1(rm_dup)', 'coverage>=4(rm_dup)', 'coverage>=10(rm_dup)', 'coverage>=20(rm_dup)']
            tmp = pd.read_csv(gene_qc, sep='\t', usecols=cols)
            tmp['sample'] = sample
            gene = os.path.basename(gene_qc).split('.')[1]
            tmp.rename(columns={col: f'{col}_{gene}' for col in cols}, inplace=True)
            if tmp_.empty:
                tmp_ = tmp.copy()
            else:
                tmp_ = pd.merge(tmp_, tmp, on=['sample'])
        df_gene_qc = df_gene_qc.append(tmp_, sort=False)
    # merge ado
    ado_dic = defaultdict(dict)
    for ado in glob.glob(f'{work_dir}/Hap/*.ADO.tsv'):
        tmp = pd.read_csv(ado, sep='\t', index_col=0)
        for row in tmp.itertuples():
            sample = row.Index
            ado_dic[sample]['sample'] = sample
            try:
                rate = round(int(row.ADO) / int(row.meet), 4)
            except (ZeroDivisionError, TypeError, ValueError):
                rate = '.'
            ado_dic[sample][f'ADO_{row.gene}'] = rate
    df_ado = pd.DataFrame.from_dict(ado_dic, orient='index')
    try:
        df = df_wgs_qc.merge(df_gene_qc.merge(df_ado, on=['sample']), on=['sample'])
    except KeyError:
        df = df_ado
    return df


def main():
    parser = ArgumentParser()
    parser.add_argument('-work_dir', help='work dir')
    parser.add_argument('-info', help='sample info')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    info = pd.read_csv(parsed_args.info, sep='\t')
    info.rename(columns={'gender': 'gender_provided'}, inplace=True)
    info.drop_duplicates(subset=['sample'], inplace=True)
    samples = info.loc[info['is_lfr'] != 'yes']['sample'].values
    qc = merge_qc(parsed_args.work_dir, samples)
    qc = pd.merge(info[['sample', 'sampleID', 'gender_provided']], qc, on=['sample'])
    qc.to_excel(parsed_args.out, index=False)


if __name__ == '__main__':
    main()
