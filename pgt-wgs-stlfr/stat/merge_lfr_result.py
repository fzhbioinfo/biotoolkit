# -*- coding:utf-8 -*-
from argparse import ArgumentParser
import pandas as pd
import glob
import os


def final_qc(lfr_dir):
    # merge longranger summary.csv
    df_summary = pd.DataFrame()
    for summary in glob.glob(f'{lfr_dir}/*/outs/summary.csv'):
        tmp = pd.read_csv(summary, sep=',')
        sample = os.path.basename(os.path.dirname(os.path.dirname(summary)))
        tmp['q30_bases_fract'] = (tmp['r1_q30_bases_fract'] + tmp['r2_q30_bases_fract']) / 2
        tmp['sample'] = sample
        cols = ['sample', 'number_reads', 'q30_bases_fract', 'bc_on_whitelist', 'mapped_reads', 'pcr_duplication',
                'snps_phased', 'n50_phase_block', 'mean_depth']
        df_summary = df_summary.append(tmp[cols])
    # merge wgs qc
    df_wgs_qc = pd.DataFrame()
    for wgs_qc in glob.glob(f'{lfr_dir}/*/Stat/*.wgs.qc.tsv'):
        tmp = pd.read_csv(wgs_qc, sep='\t')
        sample = os.path.basename(os.path.dirname(os.path.dirname(wgs_qc)))
        tmp['sample'] = sample
        cols = ['sample', 'mean_depth(rm_dup)', 'coverage>=1(rm_dup)', 'coverage>=4(rm_dup)', 'coverage>=10(rm_dup)',
                'coverage>=20(rm_dup)', 'depX', 'depY', 'XRatio', 'YRatio', 'gender']
        df_wgs_qc = df_wgs_qc.append(tmp[cols])
    # merge gene qc
    df_gene_qc = pd.DataFrame()
    for sample in ['Father', 'Mother']:
        tmp_ = pd.DataFrame()
        for gene_qc in glob.glob(f'{lfr_dir}/{sample}/Stat/*.gene.qc.tsv'):
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
    return df_summary.merge(df_wgs_qc.merge(df_gene_qc, on=['sample']), on=['sample'])


def merge_result(result_files):
    df = pd.DataFrame()
    for result in glob.glob(result_files):
        tmp = pd.read_csv(result, sep='\t')
        df = df.append(tmp, sort=False)
    return df


def main():
    parser = ArgumentParser()
    parser.add_argument('-info', help='sample info')
    parser.add_argument('-lfr_dir', help='lfr dir')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    info = pd.read_csv(parsed_args.info, sep='\t', usecols=['sample', 'sampleID', 'gender'])
    info.rename(columns={'gender': 'gender_provided'}, inplace=True)
    info.drop_duplicates(subset=['sample'], inplace=True)
    qc = final_qc(parsed_args.lfr_dir)
    ps_stat = merge_result(f'{parsed_args.lfr_dir}/*/Stat/*.PS_stat.tsv')
    ps_detect = merge_result(f'{parsed_args.lfr_dir}/*/Stat/*.PS_detect.tsv')
    qc = pd.merge(info, qc, on=['sample'])
    ps_stat = pd.merge(info[['sample', 'sampleID']], ps_stat, on=['sample'])
    ps_detect = pd.merge(info[['sample', 'sampleID']], ps_detect, on=['sample'])
    writer = pd.ExcelWriter(parsed_args.out)
    qc.to_excel(writer, index=False, sheet_name='QC')
    ps_stat.to_excel(writer, index=False, sheet_name='PS_stat')
    ps_detect.to_excel(writer, index=False, sheet_name='PS_detect')
    writer.save()
    writer.close()


if __name__ == '__main__':
    main()
