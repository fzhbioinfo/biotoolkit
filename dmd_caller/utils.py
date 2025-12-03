"""DMD caller
"""
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from functools import partial
import pandas as pd
import yaml
import sys
import os
from .statistics import cal_depth_coverage, lowess_fit, cal_cv, plot_exon


def is_continuously(df):
    """判断cnv是否连续
    """
    df = df.copy()
    df.reset_index(inplace=True, drop=True)
    df['continuously'] = 'No'
    for i in range(df.shape[0]):
        if i != df.shape[0] - 1:
            cnv_type = df.loc[i, 'type']
            cnv_type_next = df.loc[i + 1, 'type']
            if cnv_type in ['DEL', 'DUP'] and cnv_type == cnv_type_next:
                df.loc[i, 'continuously'] = 'Yes'
                df.loc[i + 1, 'continuously'] = 'Yes'
    return df


def is_nearby_threshold(conf_dic, ratio, gender):
    """判断ratio是临近阈值
    """
    if gender == 'M':
        del_abs = abs(ratio - conf_dic['del_male_x'])
        dup_abs = abs(ratio - conf_dic['dup_male_x'])
    else:
        del_abs = abs(ratio - conf_dic['del'])
        dup_abs = abs(ratio - conf_dic['dup'])
    if del_abs < conf_dic['nearby_distance'] or dup_abs < conf_dic['nearby_distance']:
        return True
    else:
        return False


def fix_continuously(df, column, conf_dic):
    """连接不连续的CNV
    """
    gender = df['gender'].values[0]
    for i in range(df.shape[0]):
        if i < df.shape[0] - 2:
            cnv_type = df.loc[i, 'type']
            cnv_type_next = df.loc[i + 1, 'type']
            ratio = df.loc[i + 1, column]
            cnv_type_next_next = df.loc[i + 2, 'type']
            if cnv_type in ['DEL', 'DUP'] and cnv_type == cnv_type_next_next and cnv_type_next == 'NORMAL' and is_nearby_threshold(conf_dic, ratio, gender):
                df.loc[i + 1, 'type'] = cnv_type
                df.loc[i, 'continuously'] = 'Yes'
                df.loc[i + 1, 'continuously'] = 'Yes'
                df.loc[i + 2, 'continuously'] = 'Yes'
    return df


def recall_nearby_threshold_cnv(df, column, conf_dic):
    """召回临近阈值CNV
    """
    cv = df['cv'].values[0]
    shape = df.shape[0]
    if cv >= conf_dic['nice_cv'] or 'DUP' not in df['type'].values or 'M' in df['gender'].values:
        return df
    else:
        df_dup = df[df['type'] == 'DUP']
        for index in df_dup.index:
            i, j = 1, 1
            # 向前搜索
            while (index - i) >= 0:
                ratio = df.loc[index - i, column]
                if df.loc[index - i, 'type'] == 'NORMAL' and abs(conf_dic['dup'] - ratio) < conf_dic['nearby_recall_distance']:
                    df.loc[index - i, 'type'] = 'DUP'
                    df.loc[index - i, 'continuously'] = 'Yes'
                    df.loc[index, 'continuously'] = 'Yes'
                else:
                    break
                i += 1
            # 向后搜索
            while (index + j) < shape:
                ratio = df.loc[index + j, column]
                if df.loc[index + j, 'type'] == 'NORMAL' and abs(conf_dic['dup'] - ratio) < conf_dic['nearby_recall_distance']:
                    df.loc[index + j, 'type'] = 'DUP'
                    df.loc[index + j,'continuously'] = 'Yes'
                    df.loc[index, 'continuously'] = 'Yes'
                else:
                    break
                j += 1
    return df


def fix_copy_number(df, column, conf_dic):
    """单外显子拷贝数修正
    """
    del_single = conf_dic['del_single']
    dup_single = conf_dic['dup_single']
    chrom = df['chrom'].values[0]
    gender = df['gender'].values[0]
    for i in range(df.shape[0]):
        continuously = df.loc[i, 'continuously']
        cnv_type = df.loc[i, 'type']
        ratio = df.loc[i, column]
        if continuously == 'Yes' or cnv_type == 'NORMAL' or (chrom == 'chrX' and gender == 'M'):
            continue
        if cnv_type == 'DEL':
            if ratio > del_single:
                df.loc[i, 'copy_num'] = 2
                df.loc[i, 'type'] = 'NORMAL'
        elif cnv_type == 'DUP':
            if ratio < dup_single:
                df.loc[i, 'copy_num'] = 2
                df.loc[i, 'type'] = 'NORMAL'
    return df


def fix_result(conf_dic, df):
    """结果修正
    """
    df = is_continuously(df)
    df = fix_continuously(df, 'ratio_batch_fix', conf_dic)
    df = recall_nearby_threshold_cnv(df, 'ratio_batch_fix', conf_dic)
    df = fix_copy_number(df, 'ratio_batch_fix', conf_dic)
    return df


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
    result = defaultdict(dict)
    for sample in bam_dic:
        # 多进程统计区域深度覆盖度
        cal_depth_coverage_func = partial(cal_depth_coverage, bam_dic[sample])
        records = [row._asdict() for row in bed.itertuples(False)]
        with Pool(cpu) as pool:
            results = pool.map(cal_depth_coverage_func, records)
        df = pd.DataFrame(results)
        # 统计结果按转录本cds的编号升序排序
        df.sort_values(by=['chrom', 'transcript', 'number'], ascending=True, inplace=True)
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
        df.index = df['chrom'] + '_' + df['transcript'] + '_' + df['exon']
        result[sample]['gc_fix'] = df
        df_gc_fix = df[['ratio_gc_fix']].copy()
        df_gc_fix.rename(columns={'ratio_gc_fix': sample}, inplace=True)
        result[sample]['gc_fix_ratio'] = df_gc_fix
        dir_name = f'{out_dir}/{sample}/coverage'
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        df.to_csv(f'{dir_name}/{sample}.gc_fix.tsv', index=False, sep='\t')
    return result


def batch_correction(bam_dic, out_dir, result, conf_dic):
    """批次校正
    """
    # 所有样本(包括control)gc校正ratio
    gc_fix_ratio = pd.concat([result[sample]['gc_fix_ratio'] for sample in result.keys()], axis=1)
    # 男女样本及其ratio参考线
    male = [s for s in result.keys() if result[s]['gender'] == 'M']
    female = [s for s in result.keys() if result[s]['gender'] == 'F']
    line_male = gc_fix_ratio[male].median(axis=1)
    line_female = gc_fix_ratio[female].median(axis=1)
    # 所有样本批次校正ratio
    line_all = gc_fix_ratio.median(axis=1)
    # control都是0, 置为1, 避免出现inf
    line_male.replace(0, 1, inplace=True)
    line_female.replace(0, 1, inplace=True)
    line_all.replace(0, 1, inplace=True)
    # 批次修正
    for sample in bam_dic:
        df = result[sample]['gc_fix']
        # 对照中同性别10个以上，该性别的样本校正用同性别的样本作为对照，否则使用全部的作为对照
        if len(male) >= 10 and result[sample]['gender'] == 'M':
            line = line_male
        elif len(female) >= 10 and result[sample]['gender'] == 'F':
            line = line_female
        else:
            line = line_all
        # 设置index与control ratio index一致
        df['ratio_batch_fix'] = df['ratio_gc_fix'] / line
        df['ratio_control'] = line
        # 拷贝数转换
        df = ratio_to_copy_number(df, 'ratio_batch_fix', conf_dic)
        # DMD波动统计
        dmd = df.loc[df['gene_symbol'] == 'DMD']
        cv = cal_cv(dmd, 'ratio_batch_fix', conf_dic['del'], conf_dic['dup'])
        df['cv'] = cv
        dir_name = f'{out_dir}/{sample}/coverage'
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        df.to_csv(f'{dir_name}/{sample}.batch_fix.tsv', index=False, sep='\t')
        result[sample]['batch_fix'] = df
    return result


def run(args):
    """dmd caller流程
    """
    # 多进程数目
    cpu = min(cpu_count(), args.cpu)
    # 阈值配置
    with open(args.conf, 'r', encoding='utf-8') as f:
        conf_dic = yaml.load(f, Loader=yaml.FullLoader)
    # 统计模式, 统计完退出
    if args.mode == 'stat':
        # 基因转录本坐标信息
        bed = pd.read_csv(args.bed, sep='\t', dtype={'entrez_id': str})
        # cds编号
        bed['number'] = bed['cds'].str.extract(r'(\d+)').astype(int)
        # cds长度
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
            df = pd.read_csv(file_path, sep='\t', dtype={'entrez_id': str})
            df.index = df['chrom'] + '_' + df['transcript'] + '_' + df['exon']
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
            df.index = df['chrom'] + '_' + df['transcript'] + '_' + df['exon']
            sample = df['sample'].values[0]
            gender = df['gender'].values[0]
            result[sample]['gender'] = gender
            df_gc_fix = df[['ratio_gc_fix']].copy()
            df_gc_fix.rename(columns={'ratio_gc_fix': sample}, inplace=True)
            result[sample]['gc_fix_ratio'] = df_gc_fix
    # 批次校正
    result = batch_correction(bam_dic, args.out_dir, result, conf_dic)
    # 画图的基因列表
    genes = pd.read_csv(args.gene, sep='\t', header=None, dtype=str)[0].tolist()
    # 修正结果函数
    fix_result_func = partial(fix_result, conf_dic)
    for sample in bam_dic:
        df = result[sample]['batch_fix']
        # 结果修正
        df_list = [df.loc[df['gene_symbol'] == gene] for gene in df['gene_symbol'].unique()]
        with Pool(cpu) as pool:
            fix_list = pool.map(fix_result_func, df_list)
        df_fix = pd.concat(fix_list, axis=0)
        df_fix.sort_values(by=['chrom', 'transcript', 'number'], ascending=True, inplace=True)
        dir_name = f'{args.out_dir}/{sample}/coverage'
        df_fix.to_csv(f'{dir_name}/{sample}.result.tsv', index=False, sep='\t')
        # 画图
        graph_dir = f'{args.out_dir}/{sample}/Exon_graph'
        if not os.path.exists(graph_dir):
            os.makedirs(graph_dir)
        for gene in genes:
            df_gene = df.loc[df['entrez_id'] == gene]
            plot_exon(df_gene, 'ratio_batch_fix', graph_dir)
