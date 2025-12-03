from argparse import ArgumentParser
from interval import Interval
from io import StringIO
import pandas as pd
import subprocess
import logging
import sys
import os


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)


# 未指定绝对路径,则环境中需要bcftools
bcftools = 'bcftools'


def get_region_ps(lfr, chromosome, start, stop):
    # 获取基因内的PS,过滤掉未分型的
    cmd = f"{bcftools} query -r {chromosome}:{start}-{stop} -f '%CHROM\t%POS\t[%PS\t%GT]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['chromosome', 'pos', 'ps', 'gt'], dtype={'ps': str, 'gt': str})
    block = block.loc[block['gt'].str.contains('\|') & (block['ps'] != '.')]
    return block['ps'].unique()


def parse_lfr(lfr):
    # lfr vcf解析,过滤掉未分型的
    cmd = f"{bcftools} query -f '%CHROM\t%POS\t%REF\t%ALT\t[%PS\t%GT]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['chromosome', 'pos', 'ref', 'alt', 'ps', 'gt'], dtype={'ps': str, 'gt': str})
    return block.loc[block['gt'].str.contains('\|') & (block['ps'] != '.')]


def get_ps_range(block, ps):
    # 获取PS范围,只保留杂合snp的
    ps_block = block.loc[block['ps'] == ps]
    return ps_block['pos'].values[0], ps_block['pos'].values[-1], ps_block.loc[(ps_block['ref'].map(len) == 1) & (ps_block['alt'].map(len) == 1) & (ps_block['gt'] != '1|1')].shape[0]


def cal_overlap_rate(ps_start, ps_stop, start, stop):
    # 求交集
    ps_interval = Interval(ps_start, ps_stop)
    gene_interval = Interval(start, stop)
    overlap = ps_interval & gene_interval
    overlap_rate = (overlap.upper_bound - overlap.lower_bound) / (stop - start)
    return overlap_rate


def main():
    parser = ArgumentParser()
    parser.add_argument('-region', help='gene region file')
    parser.add_argument('-lfr', help='lfr phased_variants.vcf.gz')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    region = pd.read_csv(parsed_args.region, sep='\t', header=None, names=['chromosome', 'start', 'stop', 'gene'])
    stat = list()
    sample = os.path.basename(os.path.dirname(os.path.dirname(parsed_args.lfr)))
    block = parse_lfr(parsed_args.lfr)
    for row in region.itertuples(index=False):
        dic = row._asdict()
        dic['sample'] = sample
        pss = get_region_ps(parsed_args.lfr, row.chromosome, row.start, row.stop)
        if len(pss) > 0:
            for ps in pss:
                dic['PS'] = ps
                dic['PS_start'], dic['PS_stop'], dic['PS_size'] = get_ps_range(block, ps)
                dic['Overlap_Rate'] = cal_overlap_rate(dic['PS_start'], dic['PS_stop'], row.start, row.stop)
        else:
            dic['PS'], dic['PS_start'], dic['PS_stop'], dic['PS_size'], dic['Overlap_Rate'] = '.', '.', '.', 0, 0
        stat.append(dic)
    stat_df = pd.DataFrame(stat)
    stat_df.loc[stat_df['PS_size'] > 1].to_csv(parsed_args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
