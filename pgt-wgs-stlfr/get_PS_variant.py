# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
from io import StringIO
import pandas as pd
import subprocess
import logging
import vcf
import sys
import os


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)


# 未指定绝对路径,则环境中需要bcftools
bcftools = 'bcftools'


def get_ps_region(lfr, chromosome, ps):
    # 提取PS的起始和终止
    cmd = f"{bcftools} query -r {chromosome} -f '%CHROM\t%POS\t[%PS]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'PS'], dtype=str)
    block = block.loc[block['PS'] == str(ps)]
    pos = block['POS'].tolist()
    return int(pos[0]), int(pos[-1])


def parse_ps_block(lfr, chromosome, start, stop):
    # bcftools解析LFR的变异分型结果,获取目标PS的block和坐标范围,并过滤无效变异
    cmd = f"{bcftools} query -r {chromosome}:{start}-{stop} -f '%CHROM\t%POS\t%REF\t[%PS\t%GT\t%TGT\t%DP\t%AD]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'PS', 'GT_lfr', 'TGT_lfr', 'DP_lfr', 'AD_lfr'],
                        dtype={'PS': str, 'DP_lfr': str})
    block = block[block.apply(lambda x: x['DP_lfr'] != '.' and int(x['DP_lfr']) > 0, axis=1)].copy()
    block['DP_lfr'] = block['DP_lfr'].astype(int)
    block['Ratio_lfr'] = (block['AD_lfr'].str.split(',').str[-1].astype(int) / block['DP_lfr']).round(3)
    return block.loc[block['GT_lfr'].str.contains('\|') & (block['PS'] != '.')]


def parse_family_block(family_vcf, chromosome, start, stop):
    family = ['Mother', 'Father']
    reader = vcf.Reader(filename=family_vcf)
    # 解析block范围内家系vcf的变异信息
    reader_ = reader.fetch(chromosome, start - 1, stop)
    stat = list()
    for record in reader_:
        # 跳过不确定的和大于2个alt的变异或低质量变异
        try:
            if 'LowQual' in record.FILTER:
                continue
        except TypeError:
            pass
        if '*' in record.ALT or len(record.ALT) > 2:
            continue
        # 每条变异信息
        call = defaultdict(dict)
        call['CHROM'] = record.CHROM
        call['POS'] = record.POS
        call['REF'] = record.REF
        alts = [str(alt) for alt in record.ALT]
        call['ALT'] = ','.join(alts)
        base = [record.REF]
        base.extend(alts)
        call['VarType'] = record.var_type
        for member in family:
            if record.genotype(member)['GT'] in ['./.', '.|.']:
                call = defaultdict(dict)
                break
            try:
                ad = record.genotype(member)['AD']
                gt = record.genotype(member)['GT'].replace('|', '/')
                dp = sum(ad)
                ratio = round(ad[-1] / dp, 3)
                call[f'{member}_GT'] = gt
                call[f'{member}_TGT'] = '/'.join([base[int(i)] for i in gt.split('/')])
                call[f'{member}_DP'] = dp
                call[f'{member}_AD'] = ','.join([str(i) for i in ad])
                call[f'{member}_Ratio'] = ratio
            except (AttributeError, ZeroDivisionError, ValueError):
                call = defaultdict(dict)
                break
        if call:
            stat.append(call)
    stat_df = pd.DataFrame(stat)
    return stat_df


def filter_aimed_variant(df, parent, hap, start, stop):
    # 过滤出可用的目标变异
    # 删除PS为空的
    df = df.loc[df['PS'].notna()]
    # 过滤snp
    df = df.loc[df['VarType'] == 'snp']
    # 过滤目标基因上下游2M内的变异
    extend = 2000000
    df = df.loc[(df['POS'] <= stop + extend) & (df['POS'] >= start - extend)].copy()
    # 按到基因中点距离对变异进行升序排序
    df['abs_distance_to_gene'] = (df['POS'] - (start + stop) / 2).abs()
    df.sort_values(by=['abs_distance_to_gene'], ascending=True, inplace=True)
    # 目标cnv是在hap1则GT_lfr筛选1|0,在hap2则筛选0|1; 
    # 看lfr是Mother/Father,那么她/他的基因型*_GT筛选0/1(与lfr同时检测为杂合更可信),另一方的基因型*_GT筛选0/0
    if parent == 'Father':
        df = df.loc[(df['Mother_GT'] == '0/0') & (df['Father_GT'] == '0/1')]
    else:
        df = df.loc[(df['Father_GT'] == '0/0') & (df['Mother_GT'] == '0/1')]
    if hap == 'hap1':
        df = df.loc[df['GT_lfr'] == '1|0']
    else:
        df = df.loc[df['GT_lfr'] == '0|1']
    # 深度均保留15X以上的
    depth = 15
    df = df.loc[(df['Mother_DP'] >= depth) & (df['Father_DP'] >= depth) & (df['DP_lfr'] >= depth)]
    # lfr检出的杂合ratio在0.4~0.6之间的
    low, up = 0.4, 0.6
    df = df.loc[(df['Ratio_lfr'] <= up) & (df['Ratio_lfr'] >= low)]
    return df


def main():
    parser = ArgumentParser()
    parser.add_argument('-family_vcf', help='family vcf file')
    parser.add_argument('-father_lfr', help='lfr phased_variants.vcf.gz')
    parser.add_argument('-mother_lfr', help='lfr phased_variants.vcf.gz')
    parser.add_argument('-hap', help='hap1 or hap2')
    parser.add_argument('-ps', help='PS number')
    parser.add_argument('-out', help='out file')
    parser.add_argument('-gene', help='gene name')
    parser.add_argument('-target', help='target.gene')
    parsed_args = parser.parse_args()
    # 获取基因范围
    target = pd.read_csv(parsed_args.target, sep='\t', header=None, usecols=[0, 1, 4, 5], names=['gene', 'chromosome', 'start', 'stop'])
    target.index = target['gene'].tolist()
    target_dic = target.to_dict('index')
    chromosome = target_dic[parsed_args.gene]['chromosome']
    gene_start = target_dic[parsed_args.gene]['start']
    gene_stop = target_dic[parsed_args.gene]['stop']
    # 确实lfr来源Father/Mother
    if parsed_args.hap not in ['hap1', 'hap2']:
        raise ValueError('hap must be hap1 or hap2!')
    lfr, parent = (parsed_args.father_lfr, 'Father') if parsed_args.father_lfr else (parsed_args.mother_lfr, 'Mother')
    # PS范围
    start, stop = get_ps_region(lfr, chromosome, parsed_args.ps)
    # PS范围家系变异
    family_block = parse_family_block(parsed_args.family_vcf, chromosome, start, stop)
    # PS范围lfr变异
    lfr_block = parse_ps_block(lfr, chromosome, start, stop)
    # 合并过滤
    family_block_lfr = pd.merge(family_block, lfr_block, on=['CHROM', 'POS', 'REF'], how='left')
    family_block_lfr_ = filter_aimed_variant(family_block_lfr, parent, parsed_args.hap, gene_start, gene_stop)
    if os.path.splitext(parsed_args.out)[-1] in ['.xls', '.xlsx']:
        family_block_lfr.to_excel(parsed_args.out, index=False)
        family_block_lfr_.to_excel(f'{parsed_args.out}.filter.xlsx', index=False)
    else:
        family_block_lfr.to_csv(parsed_args.out, index=False, sep='\t')
        family_block_lfr_.to_csv(f'{parsed_args.out}.filter.tsv', index=False, sep='\t')


if __name__ == '__main__':
    main()
