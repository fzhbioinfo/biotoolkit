# -*- coding:utf-8 -*-
from multiprocessing import cpu_count, Pool
from argparse import ArgumentParser
from pysam import AlignmentFile
from functools import reduce
import pandas as pd
import numpy as np
import logging
import os


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


ROOT = os.path.abspath(os.path.dirname(__file__))
CONFIG = os.path.join(ROOT, '../etc')


def cal_coverage(cover):
    cutoff = [1, 4, 10, 20]
    coverage = (np.sum(np.sum(cover, axis=0) >= i) for i in cutoff)
    return coverage


def statistics(args):
    bam = AlignmentFile(args[0], 'rb')
    bed = args[1]
    stat = list()
    for row in bed.itertuples(index=False):
        dic = dict()
        logger.info(f'count coverage: {row.chromosome}:{row.start}-{row.stop}')
        dic['chromosome'] = row.chromosome
        dic['start'] = row.start -1 if row.start == row.stop else row.start
        dic['stop'] = row.stop
        dic['size'] = dic['stop'] - dic['start']
        dic['reads_on_target(rm_dup)'] = bam.count(dic['chromosome'], dic['start'], dic['stop'], read_callback='all')
        dic['reads_on_target(dup)'] = bam.count(dic['chromosome'], dic['start'], dic['stop'], read_callback=lambda x: x.is_duplicate)
        cover = bam.count_coverage(dic['chromosome'], dic['start'], dic['stop'], read_callback='all')
        dic['base_on_target(rm_dup)'] = np.sum(cover)
        dic['base_on_target(dup)'] = np.sum(bam.count_coverage(dic['chromosome'], dic['start'], dic['stop'], read_callback=lambda x: x.is_duplicate))
        dic['coverage>=1(rm_dup)'], dic['coverage>=4(rm_dup)'], dic['coverage>=10(rm_dup)'], dic['coverage>=20(rm_dup)'] = cal_coverage(cover)
        stat.append(dic)
    bam.close()
    stat_df = pd.DataFrame(stat)
    return stat_df


def analysis(args):
    if args.bed:
        logger.info(f'Use merged bed: {args.bed}')
        bed = pd.read_csv(args.bed, sep='\t', header=None)
        bed.columns = ['chromosome', 'start', 'stop']
        chromosomes = bed['chromosome'].unique().tolist()
        bed_list = [bed.loc[bed['chromosome'] == chromosome] for chromosome in chromosomes]
    else:
        logger.info(f'Use split bed, bed_dir is: {args.bed_dir}, use chromosomes in chr.list')
        chromosomes = pd.read_csv(f'{CONFIG}/chr.list', sep='\t', header=None)
        chromosomes = chromosomes[0].tolist()
        bed_list = list()
        for chromosome in chromosomes:
            bed = pd.read_csv(os.path.join(args.bed_dir, f'{chromosome}.bed'), sep='\t', header=None)
            bed.columns = ['chromosome', 'start', 'stop']
            bed_list.append(bed)
    if args.bam:
        logger.info(f'Use merged bam: {args.bam}')
        bam = [args.bam for _ in chromosomes]
    else:
        logger.info(f'Use split bam, bam_dir is: {args.bam_dir}')
        bam = [os.path.join(args.bam_dir, f'{chromosome}.bqsr.bam') for chromosome in chromosomes]
    cpu = min(args.process, cpu_count())
    with Pool(cpu) as pool:
        stat = pool.map(statistics, zip(bam, bed_list))
    stat_df = reduce(lambda x, y: x.append(y, sort=False), stat)
    stat_df.to_csv(args.out, sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-bam_dir', help='split bam dir')
    parser.add_argument('-bam', help='merged bam')
    parser.add_argument('-bed', help='bed file')
    parser.add_argument('-bed_dir', help='chromosome bed dir', default=CONFIG)
    parser.add_argument('-out', help='out file', required=True)
    parser.add_argument('-process', help='process number', default=1, type=int)
    parsed_args = parser.parse_args()
    analysis(parsed_args)


if __name__ == '__main__':
    main()
