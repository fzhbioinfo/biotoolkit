"""统计
"""
import pandas as pd
import numpy as np
from pysam import AlignmentFile
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def cal_depth_coverage(bam_file, dic) -> dict:
    """计算区域深度和覆盖度
    """
    bam = AlignmentFile(bam_file, 'rb')
    reads_count = bam.count_coverage(dic['chrom'], dic['start'] - 1, dic['stop'])
    pos_depth = np.sum(reads_count, axis=0)
    base = np.sum(pos_depth)
    mean = np.mean(pos_depth)
    coverage = np.sum(pos_depth > 0) / pos_depth.shape[0]
    dic['base'] = base
    dic['mean_depth'] = mean
    dic['coverage'] = coverage
    bam.close()
    return dic


def cal_cv(df, column, low=0.75, up=1.3, ratio=0.95):
    """计算cv, 分层考虑
    """
    df_ = df.loc[df[column] > low].copy()
    if df_.empty:
        std = df[column].std()
        mean = df[column].mean()
        cv = std / mean if float(mean) != 0 else 0
    else:
        df_normal = df_.loc[df_[column] < up].copy()
        if df_normal.shape[0] / df_.shape[0] >= ratio:
            std = df_normal[column].std()
            mean = df_normal[column].mean()
        else:
            std = df_[column].std()
            mean = df_[column].mean()
        cv = std / mean
    return cv


def lowess_fit(df) -> pd.DataFrame:
    """lowess回归拟合,校正
    """
    # 待拟合数据过滤
    gc_up = 0.8
    gc_low = 0.2
    min_depth = 10
    min_size = 20
    df_filter = df.loc[(df['mean_depth'] > min_depth) & (df['size'] > min_size) & (df['gc'] > gc_low) & (df['gc'] < gc_up)].copy()
    # 再取中间90%数据进行拟合
    quantile = df_filter['ratio'].quantile([0.05, 0.95])
    df_filter = df_filter.loc[(df_filter['ratio'] > quantile.loc[0.05]) & (df_filter['ratio'] < quantile.loc[0.95])].copy()
    df_filter['lowess'] = lowess(df_filter['ratio'], df_filter['gc'], return_sorted=False)
    columns = ['transcript', 'exon']
    df_filter = pd.merge(df, df_filter[columns + ['lowess']], on=columns, how='left')
    df_filter.sort_values(by=['gc'], ascending=True, inplace=True)
    df_filter['lowess'] = df_filter['lowess'].bfill()
    df_filter['lowess'] = df_filter['lowess'].ffill()
    df_filter['lowess'] = df_filter['lowess'].fillna(1)
    df_filter['ratio_gc_fix'] = df_filter['ratio'] / df_filter['lowess']
    df = pd.merge(df, df_filter[columns + ['lowess', 'ratio_gc_fix']], on=columns, how='left')
    return df


def plot_exon(df, column, out_dir):
    """画外显子深度图, 传入单基因的DataFrame
    """
    gene = df['gene_symbol'].values[0]
    transcript = df['transcript'].values[0]
    sample = df['sample'].values[0]
    gender = df['gender'].values[0]
    chrom = df['chrom'].values[0]
    title = f'{sample}.{gene}.{transcript}'
    exon = df['exon'].values
    depth = df[column].values
    control = df['ratio_control'].values
    coverage = (df['coverage'] * 100).values
    low = 0.3 if gender == 'M' and chrom == 'chrX' else 0.75
    up = 1.5 if gender == 'M' and chrom == 'chrX' else 1.3
    _max = df[column].max()
    if _max == 0:
        _max = 1
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(exon, depth, color='blue', width=0.8, label='depth')
    if _max >= low:
        ax.plot(exon, [low] * df.shape[0], c='black', lw=0.4)
    if _max >= up:
        ax.plot(exon, [up] * df.shape[0], c='black', lw=0.4)
    ax.plot(exon, control, '*', c='black', label='control')
    ax.set_ylabel('Relative Depth')
    ax.tick_params(axis='x', bottom=False, rotation=90)
    ax.set_title(title)
    ax.legend(loc='upper left')
    ax.set_ylim(0, _max * 1.2)
    _ax = ax.twinx()
    _ax.plot(exon, coverage, 'o', c='red', ls='-', lw=0.3, label='coverage')
    _ax.set_ylabel('Coverage(%)')
    _ax.set_ylim(0, 110)
    _ax.legend(loc='upper right')
    plt.savefig(f'{out_dir}/{title}_nor_batch.png')
    plt.close()
