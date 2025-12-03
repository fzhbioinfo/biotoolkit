# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from pysam import AlignmentFile
from interval import Interval
import pandas as pd
import numpy as np
import vcf
import re


def parse_cnv_coordinate(cnv):
    # 解析cnv染色体坐标
    chromosome, start, stop = re.findall('(.*?)[:-](.*?)[:-](.*?)$', cnv)[0]
    chromosome = f'chr{chromosome}' if not chromosome.startswith('chr') else chromosome
    return chromosome, int(start), int(stop)


def parse_snv_coordinate(snv):
    # 解析snp/indel染色体坐标
    chromosome, pos, ref, alt = re.findall('(.*?)[:-](.*?)[:-](.*?)[:-](.*?)$', snv)[0]
    chromosome = f'chr{chromosome}' if not chromosome.startswith('chr') else chromosome
    return chromosome, int(pos), ref, alt


def get_cnv_ps(lfr, chromosome, start, stop):
    # 从LFR的cnv分型结果检索目标cnv,若存在多个则取交集最大的cnv
    cnv = Interval(start, stop)
    if stop - start <= 30000:
        reader = vcf.Reader(filename=f'{lfr}/outs/dels.vcf.gz')
    else:
        reader = vcf.Reader(filename=f'{lfr}/outs/large_svs.vcf.gz')
    bam = AlignmentFile(f'{lfr}/outs/phased_possorted_bam.bam', 'rb')
    sample = reader.samples[0]
    records = reader.fetch(chromosome, start, stop)
    # 默认为None
    ps, start_, stop_, gt, dp = None, None, None, None, None
    max_rate = 0
    for record in records:
        try:
            if record.INFO['PS'] is None or record.INFO['SVTYPE'] not in ['DEL']:
                continue
        except KeyError:
            continue
        try:
            cnv_detected = Interval(record.POS, record.INFO['END'])
        except KeyError:
            continue
        overlap = cnv & cnv_detected
        overlap_rate = (overlap.upper_bound - overlap.lower_bound) / (stop - start)
        if overlap_rate > max_rate:
            max_rate = overlap_rate
            ps = record.INFO['PS']
            start_, stop_, gt = record.POS, record.INFO['END'], record.genotype(sample)['GT']
            dp = cal_depth(bam, chromosome, start_, stop_)
    bam.close()
    return ps, start_, stop_, gt, dp


def cal_depth(bam, chromosome, start, stop):
    return np.round(np.sum(bam.count_coverage(chromosome, start, stop)) / (stop - start), 3)


def get_snv_ps(lfr, chromosome, pos, ref, alt):
    # 读取LFR的snp/indel变异分型结果
    reader = vcf.Reader(filename=f'{lfr}/outs/phased_variants.vcf.gz')
    sample = reader.samples[0]
    records = reader.fetch(chromosome, pos - 1, pos)
    # 默认为None
    ps, start, stop, gt, dp = None, None, None, None, None
    for record in records:
        alt_lfr = [str(a) for a in record.ALT]
        if pos == record.POS and ref == record.REF and alt in alt_lfr:
            ps, start, stop, gt, dp = record.genotype(sample)['PS'], pos, pos, record.genotype(sample)['GT'], record.genotype(sample)['DP']
            break
    return ps, start, stop, gt, dp


def main():
    parser = ArgumentParser()
    parser.add_argument('-info', help='sample info')
    parser.add_argument('-lfr_workdir', help='lfr work dir')
    parser.add_argument('-sample', help='sample')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    info = pd.read_csv(parsed_args.info, sep='\t')
    info.drop_duplicates(subset=['sample'], inplace=True)
    info.index = info['sample'].tolist()
    info_dic = info.to_dict('index')
    genes = info_dic[parsed_args.sample]['gene'].split(',')
    vars_type = info_dic[parsed_args.sample]['var_type'].split(',')
    variants = info_dic[parsed_args.sample]['variant'].split(',')
    stat = list()
    for i, gene in enumerate(genes):
        var_type = vars_type[i]
        variant = variants[i]
        dic = dict()
        dic['sample'] = parsed_args.sample
        dic['Variant'] = f'{gene}_{var_type}_{variant}'
        if var_type in ['snv', 'indel']:
            chromosome, pos, ref, alt = parse_snv_coordinate(variant)
            dic['PS'], dic['start'], dic['stop'], dic['GT'], dic['DP'] = get_snv_ps(parsed_args.lfr_workdir, chromosome, pos, ref, alt)
        elif var_type in ['cnv']:
            chromosome, start, stop = parse_cnv_coordinate(variant)
            dic['PS'], dic['start'], dic['stop'], dic['GT'], dic['DP'] = get_cnv_ps(parsed_args.lfr_workdir, chromosome, start, stop)
        else:
            pass
        stat.append(dic)
    pd.DataFrame(stat).to_csv(parsed_args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
