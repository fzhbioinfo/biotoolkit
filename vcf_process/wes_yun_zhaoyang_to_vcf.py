import pandas as pd
from pyfaidx import Fasta
import sys


def wes_zhaoyang_to_vcf(fa, chrom, start, ref, alt):
    if pd.isna(ref):
        var_type = 'ins'
    elif pd.isna(alt):
        var_type = 'del'
    elif len(ref) == 1 and len(alt) == 1:
        var_type = 'snv'
    else:
        var_type = 'delins'
    if var_type == 'snv' or var_type == 'delins':
        pos = start
    elif var_type == 'ins':
        pos = start
        base = str(fa.get_seq(chrom, start, start)).upper()
        ref = base
        alt = f'{base}{alt}'
    elif var_type == 'del':
        pos = start - 1
        base = str(fa.get_seq(chrom, pos, pos)).upper()
        ref = f'{base}{ref}'
        alt = base
    return chrom, pos, ref, alt


def vcf_judge_mutype(df: pd.DataFrame) -> pd.DataFrame:
    """vcf格式判断突变类型

    :param df: 变异信息数据框包括'CHROM', 'POS', 'REF', 'ALT'
    :type df: pd.DataFrame
    :raises ValueError: 如果是未知的突变类型
    :return: 变异信息数据框新增'MuType'
    :rtype: pd.DataFrame
    """
    # 突变类型判断, 注意delins可能有两种写法, 例: TAT -> TGGC or AT -> GGC
    df['MuType'] = 'unknown'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'snv'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() > 1), 'MuType'] = 'ins'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'del'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] == df['ALT'].str[0]), 'MuType'] = 'delins_eq'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins_neq'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins_neq'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() == 1) & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins_neq'
    if 'unknown' in df['MuType'].unique():
        raise ValueError('Unexpected MuType!')
    return df


def vcf_to_vep_format(df: pd.DataFrame) -> pd.DataFrame:
    """vcf格式转换为vep格式

    :param df: 变异信息数据框包括'CHROM', 'POS', 'REF', 'ALT', 'MuType'
    :type df: pd.DataFrame
    :return: 变异信息数据框新增'Uploaded_variation'
    :rtype: pd.DataFrame
    """
    # 根据vcf格式生成vep格式, snv与delins_neq, delins_eq一致
    df['Uploaded_variation'] = df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['REF'] + '/' + df['ALT']
    df.loc[df['MuType'] == 'ins', 'Uploaded_variation'] = df.loc[df['MuType'] == 'ins', 'CHROM'].astype(str) + '_' + (df.loc[df['MuType'] == 'ins', 'POS'] + 1).astype(str) + '_' + '-' + '/' + df.loc[df['MuType'] == 'ins', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'del', 'Uploaded_variation'] = df.loc[df['MuType'] == 'del', 'CHROM'].astype(str) + '_' + (df.loc[df['MuType'] == 'del', 'POS'] + 1).astype(str) + '_' + df.loc[df['MuType'] == 'del', 'REF'].str[1:] + '/' + '-'
    return df


def vcf_to_bgianno_format(df: pd.DataFrame) -> pd.DataFrame:
    """vcf格式转换为bgianno格式

    :param df: 变异信息数据框包括'CHROM', 'POS', 'REF', 'ALT', 'MuType'
    :type df: pd.DataFrame
    :return: 变异信息数据框新增'#Chr', 'Start', 'Stop', 'Ref', 'Call'
    :rtype: pd.DataFrame
    """
    # 根据vcf格式生成bgianno格式
    df['#Chr'] = df['CHROM']
    df['Start'] = df['POS'] - 1
    df['Stop'] = df['POS']
    df['Ref'] = df['REF']
    df['Call'] = df['ALT']
    df.loc[df['MuType'] == 'ins', 'Ref'] = '.'
    df.loc[df['MuType'] == 'ins', 'Call'] = df.loc[df['MuType'] == 'ins', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'del', 'Call'] = '.'
    df.loc[df['MuType'] == 'del', 'Ref'] = df.loc[df['MuType'] == 'del', 'REF'].str[1:]
    df.loc[df['MuType'] == 'ins', 'Start'] = df.loc[df['MuType'] == 'ins', 'Stop']
    df.loc[df['MuType'] == 'del', 'Start'] = df.loc[df['MuType'] == 'del', 'Stop']
    df.loc[df['MuType'] == 'del', 'Stop'] = df.loc[df['MuType'] == 'del', 'Stop'] + df.loc[df['MuType'] == 'del', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Stop'] = df.loc[df['MuType'] == 'delins_eq', 'Start'] + df.loc[df['MuType'] == 'delins_eq', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Start'] = df.loc[df['MuType'] == 'delins_eq', 'POS']
    df.loc[df['MuType'] == 'delins_eq', 'Ref'] = df.loc[df['MuType'] == 'delins_eq', 'REF'].str[1:]
    df.loc[df['MuType'] == 'delins_eq', 'Call'] = df.loc[df['MuType'] == 'delins_eq', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'delins_neq', 'Stop'] = df.loc[df['MuType'] == 'delins_neq', 'Start'] + df.loc[df['MuType'] == 'delins_neq', 'Ref'].str.len()
    return df


def main():
    in_file = sys.argv[1]
    fasta = Fasta(sys.argv[2])
    out_file = sys.argv[3]
    if in_file.endswith('.xlsx') or in_file.endswith('.xls'):
        df = pd.read_excel(in_file)
    else:
        df = pd.read_csv(in_file, sep='\t')
    df[['CHROM', 'POS', 'REF', 'ALT']] = df.apply(lambda x: wes_zhaoyang_to_vcf(fasta, x['chr'], x['start'], x['ref'], x['alt']), axis=1, result_type='expand')
    df = vcf_judge_mutype(df)
    df = vcf_to_vep_format(df)
    df.to_excel(out_file, index=False)
    df.to_csv(f'{out_file}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
