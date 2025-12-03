# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from interval import Interval
from io import StringIO
import pandas as pd
import subprocess
import logging
import sys


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)


# 未指定绝对路径,则环境中需要bcftools
bcftools = 'bcftools'


def get_gene_ps(lfr, chromosome, start, stop):
    # 获取基因上的PS
    cmd = f"{bcftools} query -r {chromosome}:{start}-{stop} -f '%CHROM\t%POS\t[%PS\t%GT]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'PS', 'GT'], dtype={'PS': str, 'GT': str})
    block = block.loc[block['GT'].str.contains('\|') & (block['PS'] != '.')]
    return block['PS'].unique()


def parse_lfr(lfr, chromosome):
    # 解析lfr变异
    cmd = f"{bcftools} query -r {chromosome} -f '%CHROM\t%POS\t[%PS\t%GT]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'PS', 'GT'], dtype={'PS': str, 'GT': str})
    return block.loc[block['GT'].str.contains('\|') & (block['PS'] != '.')]


def get_ps_coordinate(block, ps_list, start, stop):
    # 获取最大交集ps的区域,若无交集则输出原始基因的区域
    region = Interval(start, stop)
    max_rate = 0
    ps_start, ps_stop = start, stop
    for ps in ps_list:
        block_ = block.loc[block['PS'] == ps]
        ps_region = Interval(block_['POS'].values[0], block_['POS'].values[-1])
        overlap = region & ps_region
        overlap_rate = (overlap.upper_bound - overlap.lower_bound) / (stop - start)
        if overlap_rate >= max_rate:
            max_rate = overlap_rate
            ps_start, ps_stop = block_['POS'].values[0], block_['POS'].values[-1]
    return ps_start, ps_stop


def main():
    parser = ArgumentParser()
    parser.add_argument('-target', help='gene region file')
    parser.add_argument('-gene', help='gene symbol')
    parser.add_argument('-lfr', help='lfr phased_variants.vcf.gz')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    target = pd.read_csv(parsed_args.target, sep='\t', header=None, names=['gene', 'chromosome', 'start_extend', 'stop_extend', 'start', 'stop'])
    target.index = target['gene'].tolist()
    chromosome = target.loc[parsed_args.gene, 'chromosome']
    start = target.loc[parsed_args.gene, 'start']
    stop = target.loc[parsed_args.gene, 'stop']
    block = parse_lfr(parsed_args.lfr, chromosome)
    ps_list = get_gene_ps(parsed_args.lfr, chromosome, start, stop)
    ps_start, ps_stop = get_ps_coordinate(block, ps_list, start, stop)
    pd.DataFrame([{'chromosome': chromosome, 'ps_start': ps_start, 'ps_stop': ps_stop}]).to_csv(parsed_args.out, sep='\t', index=False, header=None)


if __name__ == '__main__':
    main()
