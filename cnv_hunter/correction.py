"""矫正
"""
import os
from functools import partial
from collections import defaultdict
from multiprocessing import Pool
import pandas as pd
import numpy as np
from .statistics import cal_depth_coverage, lowess_fit, cal_cv


def ratio_to_copy_number(df, column, conf_dic):
    """拷贝数转换
    """
    low = conf_dic['del']
    up = conf_dic['dup']
    hom_del = conf_dic['hom_del']
    df['copy_num'] = 2
    df['type'] = 'NORMAL'
    df.loc[df[column] > up, 'copy_num'] = 3
    df.loc[df[column] > up, 'type'] = 'DUP'
    df.loc[df[column] <= low, 'copy_num'] = 1
    df.loc[df[column] <= low, 'type'] = 'DEL'
    df.loc[df[column] < hom_del, 'copy_num'] = 0
    # 男性性染色体需特殊处理
    up_x = conf_dic['dup_male_x']
    low_x = conf_dic['del_male_x']
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX'), 'copy_num'] = 1
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX'), 'type'] = 'NORMAL'
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX') & (df[column] > up_x), 'copy_num'] = 2
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX') & (df[column] > up_x), 'type'] = 'DUP'
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX') & (df[column] <= low_x), 'copy_num'] = 0
    df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX') & (df[column] <= low_x), 'type'] = 'DEL'
    return df


def gc_correction(bam_dic, gender_dic, out_dir, bed, cpu):
    """深度统计与GC修正
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    result = defaultdict(dict)
    for sample in bam_dic:
        # 多进程统计区域深度覆盖度
        cal_depth_coverage_func = partial(cal_depth_coverage, bam_dic[sample])
        records = [row._asdict() for row in bed.itertuples(False)]
        with Pool(cpu) as pool:
            results = pool.map(cal_depth_coverage_func, records)
        df = pd.DataFrame(results)
        # 统计结果按坐标排序
        df.sort_values(by=['chrom', 'start', 'stop'], ascending=True, inplace=True)
        # 样本编号, 性别, 所有区域的平均深度, 每个区域的ratio
        df['sample'] = sample
        df['gender'] = gender_dic[sample]
        mean = df['base'].sum() / df['size'].sum()
        df['ratio'] = df['mean_depth'] / mean
        # 男性性染色体ratio*2
        df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX'), 'ratio'] = df.loc[(df['gender'] == 'M') & (df['chrom'] == 'chrX'), 'ratio'] * 2
        # lowess回归与gc修正
        df = lowess_fit(df)
        # gc校正等结果存储
        result[sample]['gender'] = gender_dic[sample]
        # 设置index方便批次校正除以control ratio
        df.index = df['chrom'] + '_' + df['start'].astype(str) + '_' + df['stop'].astype(str)
        result[sample]['gc_fix'] = df
        df_gc_fix = df[['ratio_gc_fix']].copy()
        df_gc_fix.rename(columns={'ratio_gc_fix': sample}, inplace=True)
        result[sample]['gc_fix_ratio'] = df_gc_fix
        df.to_csv(f'{out_dir}/{sample}.gc_fix.tsv', index=False, sep='\t')
    return result


def batch_correction(bam_dic, out_dir, result):
    """批次校正
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # 所有样本(包括control)gc校正ratio
    gc_fix_ratio = pd.concat([result[sample]['gc_fix_ratio'] for sample in result.keys()], axis=1)
    # 批次校正ratio
    line = gc_fix_ratio.median(axis=1)
    # control都是0, 置为1, 避免出现inf
    line.replace(0, 1, inplace=True)
    # 批次修正
    for sample in bam_dic:
        df = result[sample]['gc_fix']
        # 设置index与control ratio index一致
        df['ratio_batch_fix'] = df['ratio_gc_fix'] / line
        df['ratio_control'] = line
        # 计算整体cv
        cv = cal_cv(df, 'ratio_batch_fix')
        df['cv'] = cv
        # 异常值处理
        df.loc[df['ratio_batch_fix'] >= 16, 'ratio_batch_fix'] = 16
        df.loc[df['ratio_batch_fix'] <= 0.0625, 'ratio_batch_fix'] = 0.0625
        df['log2ratio'] = df['ratio_batch_fix'].apply(np.log2)
        df.to_csv(f'{out_dir}/{sample}.batch_fix.tsv', index=False, sep='\t')
