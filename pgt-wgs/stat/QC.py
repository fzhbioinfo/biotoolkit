# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
import numpy as np
import glob
import json
import os


class QC:
    def __init__(self, coverage):
        self.coverage = pd.read_csv(coverage, sep='\t')
        self.filter_dir = f'{os.path.dirname(os.path.dirname(coverage))}/Filter'

    def statistics(self, args):
        stat = defaultdict(dict)
        # fq统计
        total_reads = 0
        raw_bases = 0
        fqs = glob.glob(f'{self.filter_dir}/*/*.fastp.json')
        for fq in fqs:
            with open(fq, 'r') as f:
                f_dict = json.load(f)
                total_reads += f_dict['summary']['after_filtering']['total_reads']
                raw_bases += f_dict['summary']['before_filtering']['total_bases']
        mapped_reads = self.coverage['reads_on_target(dup)'].sum() + self.coverage['reads_on_target(rm_dup)'].sum()
        mapped_reads_rm_dup = self.coverage['reads_on_target(rm_dup)'].sum()
        # 是否包括线粒体
        if not args.use_mt:
            self.coverage = self.coverage.loc[~self.coverage['chromosome'].str.contains('M')]
        coverage = self.coverage.copy()
        coverage.drop(columns=['chromosome', 'start', 'stop'], inplace=True)
        dic = coverage.sum().to_dict()
        stat['genome_size'] = f'{round(dic["size"] / np.power(10,6))}M'
        stat['raw_bases'] = f'{round(raw_bases / np.power(10, 9))}G'
        stat['mapped_rate'] = np.round(mapped_reads / total_reads, 4)
        stat['mapped_rate(rm_dup)'] = np.round(mapped_reads_rm_dup / total_reads, 4)
        stat['dup_rate'] = np.round(dic['reads_on_target(dup)'] / (dic['reads_on_target(dup)'] + dic['reads_on_target(rm_dup)']), 4)
        stat['mean_depth(rm_dup)'] = np.round(dic['base_on_target(rm_dup)'] / dic['size'], 4)
        stat['coverage>=1(rm_dup)'] = np.round(dic['coverage>=1(rm_dup)'] / dic['size'], 4)
        stat['coverage>=4(rm_dup)'] = np.round(dic['coverage>=4(rm_dup)'] / dic['size'], 4)
        stat['coverage>=10(rm_dup)'] = np.round(dic['coverage>=10(rm_dup)'] / dic['size'], 4)
        stat['coverage>=20(rm_dup)'] = np.round(dic['coverage>=20(rm_dup)'] / dic['size'], 4)
        return stat

    @staticmethod
    def xy_ratio(x, y):
        if y > 0.25:
            gender = 'M'
        else:
            gender = 'F'
        return gender

    def analysis_gender(self):
        # autosomal = self.coverage.loc[~(self.coverage['chromosome'].isin(['chrX', 'chrY']) | self.coverage['chromosome'].str.contains('M'))].copy()
        autosomal = self.coverage.loc[~self.coverage['chromosome'].str.contains('M')].copy()
        autosomal.drop(columns=['chromosome', 'start', 'stop'], inplace=True)
        dic = autosomal.sum().to_dict()
        autosomal_mean_depth = dic['base_on_target(rm_dup)'] / dic['size']
        sex_chromosome = self.coverage.loc[self.coverage['chromosome'].isin(['chrX', 'chrY'])].copy()
        sex_chromosome_stat = sex_chromosome[['chromosome', 'size', 'base_on_target(rm_dup)']].groupby(['chromosome']).agg(np.sum)
        x_mean_depth = np.round(sex_chromosome_stat.loc['chrX', 'base_on_target(rm_dup)'] / sex_chromosome_stat.loc['chrX', 'size'], 4)
        y_mean_depth = np.round(sex_chromosome_stat.loc['chrY', 'base_on_target(rm_dup)'] / sex_chromosome_stat.loc['chrY', 'size'], 4)
        if autosomal_mean_depth == 0:
            x_ratio = '.'
            y_ratio = '.'
            gender = '.'
        else:
            x_ratio = np.round(x_mean_depth / autosomal_mean_depth, 4)
            y_ratio = np.round(y_mean_depth / autosomal_mean_depth, 4)
            gender = self.xy_ratio(x_ratio, y_ratio)
        return x_mean_depth, y_mean_depth, x_ratio, y_ratio, gender

    @classmethod
    def analysis(cls, args):
        qc = cls(args.coverage)
        stat = qc.statistics(args)
        if args.analysis_gender:
            stat['depX'], stat['depY'], stat['XRatio'], stat['YRatio'], stat['gender'] = qc.analysis_gender()
        stat_df = pd.DataFrame.from_dict(stat, orient='index').T
        stat_df.to_csv(args.out, sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-use_mt', help='use mt or not', action='store_true')
    parser.add_argument('-coverage', help='coverage file')
    parser.add_argument('-analysis_gender', help='analysis_gender', action='store_true')
    parser.add_argument('-out', help='out file', required=True)
    parsed_args = parser.parse_args()
    QC.analysis(parsed_args)


if __name__ == '__main__':
    main()
