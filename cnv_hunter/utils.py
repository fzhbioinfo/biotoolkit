"""cnv caller
"""
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import sys
import os
import yaml
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from .statistics import slide_window_gc
from .correction import batch_correction, gc_correction, ratio_to_copy_number
from .plot import plot_segment


ROOT = os.path.dirname(os.path.abspath(__file__))
CBS = os.path.join(ROOT, 'CBS.R')
RUN_CBS = f'Rscript {CBS}'


def join_segment(ratio_file, segment_file):
    """分段结果合并
    """
    df = pd.read_csv(ratio_file, sep='\t')
    segment = pd.read_csv(segment_file, sep='\t', usecols=['chrom', 'loc.start', 'seg.mean', 'pval'])
    rename = {'loc.start': 'start'}
    segment.rename(columns=rename, inplace=True)
    df = df.merge(segment, on=['chrom', 'start'], how='left')
    df['seg.mean'] = df['seg.mean'].ffill()
    df['copy_ratio'] = df['seg.mean'].apply(np.exp2)
    return df


def run(args):
    """cnv caller流程
    """
    # 多进程数目
    cpu = min(cpu_count(), args.cpu)
    # 滑动窗模式, 统计完退出
    if args.mode == 'slide_window':
        bed = pd.read_csv(args.bed, sep='\t', header=None, names=['chrom', 'start', 'stop'])
        bed = bed.loc[bed['chrom'].str.contains('Y|M') == False]
        fasta = Fasta(args.fasta)
        slide_window_gc(bed, args.win, args.slide, args.win_min, fasta).to_csv(args.out, index=False, sep='\t')
        sys.exit(0)
    # 深度统计和gc校正
    elif args.mode == 'stat':
        # slide window模式下坐标信息
        bed = pd.read_csv(args.bed, sep='\t')
        bed['size'] = bed['stop'] - bed['start'] + 1
        # 性别信息
        gender = pd.read_csv(args.gender, sep='\t', header=None, dtype=str, names=['sample', 'gender'])
        gender_dic = gender.set_index('sample')['gender'].to_dict()
        # bam路径
        bam = pd.read_csv(args.bam, sep='\t', header=None, dtype=str, names=['sample', 'bam'])
        bam_dic = bam.set_index('sample')['bam'].to_dict()
        # 存储结果字典, 增加固定control时去更新/增加样本的gc_fix_ratio和gender
        result = gc_correction(bam_dic, gender_dic, args.out_dir, bed, cpu)
        sys.exit(0)
    else:
        # 非统计模式需要读取事先gc修正等统计结果
        result = defaultdict(dict)
        bam_dic = {}
        gender_dic = {}
        gc_fix_files = pd.read_csv(args.gc_fix, sep='\t', header=None, dtype=str)[0].tolist()
        for file_path in gc_fix_files:
            df = pd.read_csv(file_path, sep='\t')
            df.index = df['chrom'] + '_' + df['start'].astype(str) + '_' + df['stop'].astype(str)
            sample = df['sample'].values[0]
            gender = df['gender'].values[0]
            result[sample]['gender'] = gender
            result[sample]['gc_fix'] = df
            df_gc_fix = df[['ratio_gc_fix']].copy()
            df_gc_fix.rename(columns={'ratio_gc_fix': sample}, inplace=True)
            result[sample]['gc_fix_ratio'] = df_gc_fix
            bam_dic[sample] = None
            gender_dic[sample] = gender
    # 如果有固定control, 注意control样本编号与待测样本编号不能相同
    if args.control:
        control_files = pd.read_csv(args.control, sep='\t', header=None, dtype=str)[0].tolist()
        for file_path in control_files:
            df = pd.read_csv(file_path, sep='\t')
            df.index = df['chrom'] + '_' + df['start'].astype(str) + '_' + df['stop'].astype(str)
            sample = df['sample'].values[0]
            gender = df['gender'].values[0]
            df_gc_fix = df[['ratio_gc_fix']].copy()
            df_gc_fix.rename(columns={'ratio_gc_fix': sample}, inplace=True)
            result[sample]['gc_fix_ratio'] = df_gc_fix
            result[sample]['gender'] = gender
    # 批次校正
    batch_correction(bam_dic, args.out_dir, result)
    # segment
    commands = [f'{RUN_CBS} {args.out_dir}/{sample}.batch_fix.tsv {args.out_dir}/{sample}.segment.tsv' for sample in bam_dic]
    with Pool(cpu) as pool:
        pool.map(os.system, commands)
    # 配置
    with open(args.conf, 'r', encoding='utf-8') as f:
        conf_dic = yaml.load(f, Loader=yaml.FullLoader)
    genome = 'GRCh38' if '38' in args.genome else 'GRCh37'
    cnv_name = pd.read_csv(os.path.join(ROOT, conf_dic[f'cnv_name_{genome}']), sep='\t', header=None, names=['chrom', 'start', 'stop', 'name'])
    # 合并结果和画图
    for sample in bam_dic:
        df = join_segment(f'{args.out_dir}/{sample}.batch_fix.tsv', f'{args.out_dir}/{sample}.segment.tsv')
        df = ratio_to_copy_number(df, 'copy_ratio', conf_dic)
        df.to_csv(f'{args.out_dir}/{sample}.result.tsv', index=False, sep='\t')
        plot_segment(df, cnv_name, 'ratio_batch_fix', f'{args.out_dir}/{sample}.jpg')
