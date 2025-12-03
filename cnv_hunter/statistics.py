"""统计
"""
import re
import pandas as pd
import numpy as np
from pysam import AlignmentFile
from statsmodels.nonparametric.smoothers_lowess import lowess


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


def lowess_fit(df) -> pd.DataFrame:
    """lowess回归拟合,校正
    """
    # 待拟合数据过滤
    gc_up = 0.8
    gc_low = 0.2
    min_depth = 10
    min_size = 20
    max_ratio = 5
    min_ratio = 0.1
    df_filter = df.loc[(df['mean_depth'] > min_depth) & (df['size'] > min_size) & (df['gc'] > gc_low) & (df['gc'] < gc_up) & (df['ratio'] < max_ratio) & (df['ratio'] > min_ratio)].copy()
    # 再取中间90%数据进行拟合
    quantile = df_filter['ratio'].quantile([0.05, 0.95])
    df_filter = df_filter.loc[(df_filter['ratio'] > quantile.loc[0.05]) & (df_filter['ratio'] < quantile.loc[0.95])].copy()
    df_filter['lowess'] = lowess(df_filter['ratio'], df_filter['gc'], return_sorted=False)
    columns = ['chrom', 'start', 'stop']
    df_filter = pd.merge(df, df_filter[columns + ['lowess']], on=columns, how='left')
    df_filter.sort_values(by=['gc'], ascending=True, inplace=True)
    df_filter['lowess'] = df_filter['lowess'].bfill()
    df_filter['lowess'] = df_filter['lowess'].ffill()
    df_filter['lowess'] = df_filter['lowess'].fillna(1)
    df_filter['ratio_gc_fix'] = df_filter['ratio'] / df_filter['lowess']
    df = pd.merge(df, df_filter[columns + ['lowess', 'ratio_gc_fix']], on=columns, how='left')
    return df


def slide_window(start, stop, win, slide, win_min):
    """
    :param start: 起始坐标
    :param stop: 终止坐标
    :param win: 窗口大小
    :param slide: 窗口滑动大小
    :param win_min: stop - start < win_min则忽略此区域
    :return: 滑动窗口区域坐标列表（划到最后剩余的长度不足一个窗口长度时，并入最后一个窗口）
    """
    intervals = list()
    _stop = start + win
    if stop - start < win_min:
        return None
    if _stop >= stop:
        intervals.append((start, stop))
        return intervals
    while True:
        intervals.append((start, _stop))
        start += slide
        _stop = start + win
        if _stop > stop:
            intervals.pop()
            start -= slide
            intervals.append((start, stop))
            break
    return intervals


def cal_region_gc_content(fasta, chrom, start, stop):
    """
    计算区域的GC含量, N不包括在内
    """
    seq = str(fasta.get_seq(chrom, start, stop))
    seq_len = len(re.findall('[^N]', seq))
    if seq_len == 0:
        return 0
    else:
        return round(len(re.findall('[gcGC]', seq)) / seq_len, 3)


def slide_window_gc(bed, win, slide, win_min, fasta) -> pd.DataFrame:
    """滑动窗口gc统计
    """
    result = []
    for row in bed.itertuples(False):
        chrom = row.chrom
        start = row.start
        stop = row.stop
        windows = slide_window(start, stop, win, slide, win_min)
        if windows:
           for win_start, win_stop in windows:
               result.append({'chrom': chrom, 'start': win_start, 'stop': win_stop, 'gc': cal_region_gc_content(fasta, chrom, win_start, win_stop)})
    df = pd.DataFrame(result)
    df.sort_values(by=['chrom', 'start', 'stop'], inplace=True)
    return pd.DataFrame(result)


def cal_cv(df, column):
    """计算CV
    """
    low, up = 0.3, 4
    df = df.loc[(df[column] > low) & (df[column] < up)]
    std = df[column].std()
    mean = df[column].mean()
    cv = round(std / mean, 3) if mean else 0
    return cv
