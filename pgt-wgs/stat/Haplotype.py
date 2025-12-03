# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
from io import BytesIO
from PIL import Image
import pandas as pd
import numpy as np
import logging
import vcf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)


class Haplotype:
    def __init__(self, vcf_file, target, gene, proband=None):
        # 类初始化,根据vcf获取家系成员,家系模式,胚胎列表
        self.reader = vcf.Reader(filename=vcf_file)
        self.embryos = [int(sample.split('Embryo')[1]) for sample in self.reader.samples if sample.startswith('Embryo')]
        self.embryos.sort()
        self.embryos = [f'Embryo{i}' for i in self.embryos]
        if proband:
            self.family = ['Mother', 'Father', proband]
            self.mode = 'core'
            self.embryos.remove(proband)
        else:
            if 'Maternal' in self.reader.samples:
                self.family = ['Mother', 'Father', 'Maternal']
                self.mode = 'non-core-Mat'
            elif 'Paternal' in self.reader.samples:
                self.family = ['Mother', 'Father', 'Paternal']
                self.mode = 'non-core-Pat'
            else:
                self.family = ['Mother', 'Father', 'Proband']
                self.mode = 'core'
        # 获取基因信息
        self.gene = gene
        self.target = pd.read_csv(target, sep='\t', header=None)
        self.target.columns = ['gene', 'chromosome', 'start_extend', 'end_extend', 'start', 'end']
        self.target.index = self.target['gene'].tolist()
        self.chromosome = self.target.loc[gene, 'chromosome']
        self.start = self.target.loc[gene, 'start']
        self.end = self.target.loc[gene, 'end']
        self.start_extend = self.target.loc[gene, 'start_extend']
        self.end_extend = self.target.loc[gene, 'end_extend']

    def parsing(self, args):
        # 获取目标区域位点并解析
        stat = list()
        reader = self.reader.fetch(self.chromosome, self.start_extend, self.end_extend)
        for record in reader:
            # 任意家系成员基因型分析失败该变异记录舍弃
            skip = False
            # 跳过不确定的和2个alt的变异或低质量变异
            try:
                if 'LowQual' in record.FILTER:
                    continue
            except TypeError:
                pass
            if '*' in record.ALT or len(record.ALT) > 1:
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
            for member in self.family:
                if record.genotype(member)['GT'] in ['./.', '.|.']:
                    call = defaultdict(dict)
                    skip = True
                    break
                try:
                    ad = record.genotype(member)['AD']
                    gt = record.genotype(member)['GT'].replace('|', '/')
                    if gt not in ['0/0', '0/1', '1/1']:
                        logger.info(f'Unexpected condition. {record.CHROM}-{record.POS}-{record.REF}-{call["ALT"]}-{gt}-{member}')
                        raise ValueError
                    dp = sum(ad)
                    ratio = round(ad[1] / dp, 3)
                    call[f'{member}_DP'] = dp
                    call[f'{member}_AD'] = ','.join([str(i) for i in ad])
                    call[f'{member}_ratio'] = ratio
                    call[f'{member}_GT_origin'] = gt
                    call[f'{member}_origin'] = '/'.join([base[int(i)] for i in gt.split('/')])
                    call[f'{member}_GT_correct'] = self.genotype_correct(gt, ratio, args.method_trios, args.ref_het_ratio, args.het_hom_ratio)
                    call[f'{member}_correct'] = '/'.join([base[int(i)] for i in call[f'{member}_GT_correct'].split('/')])
                except (AttributeError, ZeroDivisionError, ValueError):
                    call = defaultdict(dict)
                    skip = True
                    break
            for embryo in self.embryos:
                if skip:
                    break
                try:
                    if record.genotype(embryo)['GT'] in ['./.', '.|.']:
                        raise ValueError
                    ad = record.genotype(embryo)['AD']
                    gt = record.genotype(embryo)['GT'].replace('|', '/')
                    if gt not in ['0/0', '0/1', '1/1']:
                        logger.info(f'Unexpected condition. {record.CHROM}-{record.POS}-{record.REF}-{call["ALT"]}-{gt}-{embryo}')
                        raise ValueError
                    dp = sum(ad)
                    ratio = round(ad[1] / dp, 3)
                    call[f'{embryo}_DP'] = dp
                    call[f'{embryo}_AD'] = ','.join([str(i) for i in ad])
                    call[f'{embryo}_ratio'] = ratio
                    call[f'{embryo}_GT_origin'] = gt
                    call[f'{embryo}_GT_correct'] = self.genotype_correct(gt, ratio, args.method_embryo)
                    call[f'{embryo}_origin'] = '/'.join([base[int(i)] for i in gt.split('/')])
                    call[f'{embryo}_correct'] = '/'.join([base[int(i)] for i in call[f'{embryo}_GT_correct'].split('/')])
                except (AttributeError, ZeroDivisionError, ValueError):
                    call[f'{embryo}_DP'] = '.'
                    call[f'{embryo}_AD'] = '.'
                    call[f'{embryo}_ratio'] = '.'
                    call[f'{embryo}_GT_origin'] = '.'
                    call[f'{embryo}_GT_correct'] = '.'
                    call[f'{embryo}_origin'] = '.'
                    call[f'{embryo}_correct'] = '.'
            if call:
                stat.append(call)
        stat_df = pd.DataFrame(stat)
        return stat_df

    @staticmethod
    def genotype_correct(gt, ratio, method, low=0.05, up=0.95):
        # 基因型矫正
        if method == 'simple':
            if 0 < ratio < 1:
                gt_correct = '0/1'
            else:
                gt_correct = gt
        else:
            if gt == '0/1':
                if ratio < low:
                    gt_correct = '0/0'
                elif ratio > up:
                    gt_correct = '1/1'
                else:
                    gt_correct = gt
            elif gt == '1/1':
                if ratio <= up:
                    gt_correct = '0/1'
                else:
                    gt_correct = gt
            else:
                if ratio >= low:
                    gt_correct = '0/1'
                else:
                    gt_correct = gt
        return gt_correct

    @staticmethod
    def core_family_haplotype(stat_dic, suffix, proband):
        # 核心家系夫妻单体型分析
        if stat_dic[f'Father_GT_{suffix}'] == '0/0' and stat_dic[f'Mother_GT_{suffix}'] == '0/1':
            if stat_dic[f'{proband}_GT_{suffix}'] == '0/0':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['REF'], stat_dic['REF'], stat_dic['ALT'], 'yes'
            elif stat_dic[f'{proband}_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['REF'], stat_dic['ALT'], stat_dic['REF'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        elif stat_dic[f'Father_GT_{suffix}'] == '0/1' and stat_dic[f'Mother_GT_{suffix}'] == '0/0':
            if stat_dic[f'{proband}_GT_{suffix}'] == '0/0':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['ALT'], stat_dic['REF'], stat_dic['REF'], 'yes'
            elif stat_dic[f'{proband}_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['REF'], stat_dic['REF'], stat_dic['REF'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        elif stat_dic[f'Father_GT_{suffix}'] == '0/1' and stat_dic[f'Mother_GT_{suffix}'] == '1/1':
            if stat_dic[f'{proband}_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['ALT'], stat_dic['ALT'], stat_dic['ALT'], 'yes'
            elif stat_dic[f'{proband}_GT_{suffix}'] == '1/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['REF'], stat_dic['ALT'], stat_dic['ALT'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        elif stat_dic[f'Father_GT_{suffix}'] == '1/1' and stat_dic[f'Mother_GT_{suffix}'] == '0/1':
            if stat_dic[f'{proband}_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['ALT'], stat_dic['REF'], stat_dic['ALT'], 'yes'
            elif stat_dic[f'{proband}_GT_{suffix}'] == '1/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['ALT'], stat_dic['ALT'], stat_dic['REF'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        else:
            stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        return stat_dic

    @staticmethod
    def non_core_family_haplotype(stat_dic, suffix, member):
        # 非核心家系夫妻单体型分析
        if member == 'Maternal':
            if stat_dic[f'Mother_GT_{suffix}'] == '0/1' and stat_dic[f'{member}_GT_{suffix}'] == '0/0' and stat_dic[f'Father_GT_{suffix}'] == '0/0':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['REF'], stat_dic['ALT'], stat_dic['REF'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '0/1' and stat_dic[f'{member}_GT_{suffix}'] == '0/0' and stat_dic[f'Father_GT_{suffix}'] == '1/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['ALT'], stat_dic['ALT'], stat_dic['REF'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '0/1' and stat_dic[f'{member}_GT_{suffix}'] == '1/1' and stat_dic[f'Father_GT_{suffix}'] == '0/0':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['REF'], stat_dic['REF'], stat_dic['ALT'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '0/1' and stat_dic[f'{member}_GT_{suffix}'] == '1/1' and stat_dic[f'Father_GT_{suffix}'] == '1/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['ALT'], stat_dic['REF'], stat_dic['ALT'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        else:
            if stat_dic[f'Mother_GT_{suffix}'] == '0/0' and stat_dic[f'{member}_GT_{suffix}'] == '0/0' and stat_dic[f'Father_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['REF'], stat_dic['REF'], stat_dic['REF'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '1/1' and stat_dic[f'{member}_GT_{suffix}'] == '0/0' and stat_dic[f'Father_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['ALT'], stat_dic['REF'], stat_dic['ALT'], stat_dic['ALT'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '0/0' and stat_dic[f'{member}_GT_{suffix}'] == '1/1' and stat_dic[f'Father_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['ALT'], stat_dic['REF'], stat_dic['REF'], 'yes'
            elif stat_dic[f'Mother_GT_{suffix}'] == '1/1' and stat_dic[f'{member}_GT_{suffix}'] == '1/1' and stat_dic[f'Father_GT_{suffix}'] == '0/1':
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = stat_dic['REF'], stat_dic['ALT'], stat_dic['ALT'], stat_dic['ALT'], 'yes'
            else:
                stat_dic['F1'], stat_dic['F2'], stat_dic['M1'], stat_dic['M2'], stat_dic['IF'] = '.', '.', '.', '.', 'no'
        return stat_dic

    @staticmethod
    def which_allele(dic):
        # 碱基唯一的allele
        allele = [key for key in dic.keys() if list(dic.values()).count(dic[key]) == 1]
        return allele[0], dic[allele[0]]

    def haplotyping(self, args, stat_df):
        haplotype = list()
        if args.proband is None:
            proband = 'Proband'
        else:
            proband = args.proband
        # 迭代每条变异统计的信息
        for stat in stat_df.itertuples(index=False):
            # 转为OrderDict方便使用
            stat_dic = stat._asdict()
            # 夫妻单体型构建
            if self.mode == 'core':
                if args.correct in ['both', 'trios']:
                    stat_dic = self.core_family_haplotype(stat_dic, 'correct', proband)
                else:
                    stat_dic = self.core_family_haplotype(stat_dic, 'origin', proband)
            elif self.mode == 'non-core-Mat':
                if args.correct in ['both', 'trios']:
                    stat_dic = self.non_core_family_haplotype(stat_dic, 'correct', 'Maternal')
                else:
                    stat_dic = self.non_core_family_haplotype(stat_dic, 'origin', 'Maternal')
            else:
                if args.correct in ['both', 'trios']:
                    stat_dic = self.non_core_family_haplotype(stat_dic, 'correct', 'Paternal')
                else:
                    stat_dic = self.non_core_family_haplotype(stat_dic, 'origin', 'Paternal')
            # IF能够确定是哪条allele,IF添加过滤标签:家系成员该IF深度低于阈值标签为LowDepth
            dic = {'F1': stat_dic['F1'], 'F2': stat_dic['F2'], 'M1': stat_dic['M1'], 'M2': stat_dic['M2']}
            if stat_dic['IF'] == 'yes':
                stat_dic['Allele'] = self.which_allele(dic)[0]
                stat_dic['IF_filter'] = 'Pass'
                if self.mode == 'core':
                    if np.sum(np.array([stat_dic['Father_DP'], stat_dic['Mother_DP'], stat_dic[f'{proband}_DP']]) < args.triosDepth):
                        stat_dic['IF_filter'] = 'LowDepth'
                elif self.mode == 'non-core-Mat':
                    if np.sum(np.array([stat_dic['Father_DP'], stat_dic['Mother_DP'], stat_dic['Maternal_DP']]) < args.triosDepth):
                        stat_dic['IF_filter'] = 'LowDepth'
                else:
                    if np.sum(np.array([stat_dic['Father_DP'], stat_dic['Mother_DP'], stat_dic['Paternal_DP']]) < args.triosDepth):
                        stat_dic['IF_filter'] = 'LowDepth'
            else:
                stat_dic['Allele'] = '.'
                stat_dic['IF_filter'] = 'Fail'
            # 判断变异的位置
            if stat_dic['POS'] < self.start:
                stat_dic['Location'] = 'upstream'
            elif stat_dic['POS'] > self.end:
                stat_dic['Location'] = 'downstream'
            else:
                stat_dic['Location'] = 'internal'
            # 胚胎遗传的单体型分析
            for embryo in self.embryos:
                if args.correct in ['both', 'embryo']:
                    embryo_gt = stat_dic[f'{embryo}_GT_correct']
                    embryo_ = stat_dic[f'{embryo}_correct']
                else:
                    embryo_gt = stat_dic[f'{embryo}_GT_origin']
                    embryo_ = stat_dic[f'{embryo}_origin']
                if stat_dic['IF'] == 'no' or embryo_gt == '.' or (len(set(embryo_.split('/'))) == 1 and embryo_.split('/')[0] != self.which_allele(dic)[1]):
                    stat_dic[f'{embryo}_hap'] = '.'
                elif stat_dic[f'{embryo}_DP'] < args.embryoDepth:
                    stat_dic[f'{embryo}_hap'] = f'{stat_dic["Allele"]}-LowDepth'
                else:
                    stat_dic[f'{embryo}_hap'] = f'{stat_dic["Allele"]}'
            haplotype.append(stat_dic)
        haplotype_df = pd.DataFrame(haplotype)
        if haplotype_df.empty:
            columns = ['CHROM', 'POS', 'REF', 'ALT', 'VarType'] + [f'{sample}_{col}' for sample in self.reader.samples for col in ['DP', 'AD', 'ratio', 'GT_origin', 'origin', 'GT_correct', 'correct']] + ['F1', 'F2', 'M1', 'M2', 'IF', 'Allele', 'IF_filter', 'Location'] + [f'{embryo}_hap' for embryo in self.embryos]
            haplotype_df = pd.DataFrame(columns=columns)
        haplotype_df.to_csv(f'{args.prefix}.{self.gene}.hap.tsv', index=False, sep='\t')
        return haplotype_df

    @staticmethod
    def counting(upstream, internal, downstream, sample, gene):
        # 支持的单体型计数
        count = defaultdict(dict)
        count['sample'] = sample
        count['Gene'] = gene
        if sample == 'IF':
            upstream_count = upstream['Allele'].value_counts()
            internal_count = internal['Allele'].value_counts()
            downstream_count = downstream['Allele'].value_counts()
        else:
            upstream_count = upstream[f'{sample}_hap'].value_counts()
            internal_count = internal[f'{sample}_hap'].value_counts()
            downstream_count = downstream[f'{sample}_hap'].value_counts()
        for i in ['F1', 'F2', 'M1', 'M2']:
            try:
                count[f'{i}-U'] = upstream_count.loc[i]
            except KeyError:
                count[f'{i}-U'] = 0
            try:
                count[f'{i}-I'] = internal_count.loc[i]
            except KeyError:
                count[f'{i}-I'] = 0
            try:
                count[f'{i}-D'] = downstream_count.loc[i]
            except KeyError:
                count[f'{i}-D'] = 0
        return count

    @staticmethod
    def resulting(score_df, embryos, up=0.2, low=0.05):
        # 单体型判断结果,目标单体型IF比例0.2,非目标0.05
        score_df.index = score_df['sample'].tolist()
        for embryo in embryos:
            for fm in ['F', 'M']:
                result = list()
                for uid in ['U', 'I', 'D']:
                    if_1st = score_df.loc['IF', f'{fm}1-{uid}']
                    if_2nd = score_df.loc['IF', f'{fm}2-{uid}']
                    embryo_1st = score_df.loc[embryo, f'{fm}1-{uid}']
                    embryo_2nd = score_df.loc[embryo, f'{fm}2-{uid}']
                    if if_1st > 0 and if_2nd > 0:
                        r1 = embryo_1st / if_1st
                        r2 = embryo_2nd / if_2nd
                        if r1 > up and r2 < low:
                            result.append(f'{fm}1')
                        elif r1 < low and r2 > up:
                            result.append(f'{fm}2')
                        else:
                            result.append('-')
                    elif if_1st > 0 and if_2nd == 0:
                        r1 = embryo_1st / if_1st
                        if r1 > up:
                            result.append(f'{fm}1')
                        else:
                            result.append('-')
                    elif if_1st == 0 and if_2nd > 0:
                        r2 = embryo_2nd / if_2nd
                        if r2 > up:
                            result.append(f'{fm}2')
                        else:
                            result.append('-')
                    else:
                        result.append('-')
                if fm == 'F':
                    score_df.loc[embryo, 'Paternal'] = '|'.join(result)
                else:
                    score_df.loc[embryo, 'Maternal'] = '|'.join(result)
        return score_df

    def scoring(self, args, haplotype_df):
        # IF过滤单体型分型结果输出到表格和结果可视化
        score = list()
        if not args.keep_indel:
            haplotype_df = haplotype_df[haplotype_df['VarType'] != 'indel']
        haplotype_df = haplotype_df[haplotype_df['IF_filter'] == 'Pass']
        upstream = haplotype_df[haplotype_df['POS'] < self.start]
        internal = haplotype_df[(haplotype_df['POS'] >= self.start) & (haplotype_df['POS'] <= self.end)]
        downstream = haplotype_df[haplotype_df['POS'] > self.end]
        for sample in ['IF'] + self.embryos:
            score.append(self.counting(upstream, internal, downstream, sample, self.gene))
        score_df = pd.DataFrame(score)
        score_df['Paternal'], score_df['Maternal'] = '', ''
        if self.embryos:
            score_df = self.resulting(score_df, self.embryos)
            self.scatter(args, haplotype_df, self.gene)
        score_df.to_excel(f'{args.prefix}.{self.gene}.score.xlsx', index=False)

    def scatter(self, args, haplotype_df, gene):
        # 单体型结果画图
        imgs = list()
        ticks = ['F1', 'F2', 'M1', 'M2']
        y_locator = plt.MultipleLocator(1)
        for embryo in self.embryos:
            fig, ax = plt.subplots(figsize=(10, 3))
            ax.yaxis.set_major_locator(y_locator)
            ax.set_ylim(0.5, 4.5)
            plt.yticks(range(1, 5), ticks)
            ax.set_title(embryo, {'fontweight': 'bold'})
            hap = haplotype_df[haplotype_df[f'{embryo}_hap'].isin(ticks)]
            start = hap.loc[hap['POS'] < self.start].shape[0]
            end = hap.loc[hap['POS'] < self.end].shape[0]
            ax.axvline(x=start - 1, ls=':', c='black')
            ax.axvline(x=end - 1, ls=':', c='black')
            for i, j in enumerate(ticks):
                values = (hap[f'{embryo}_hap'] == j).replace(True, i+1)
                values = values.replace(0, np.nan)
                ax.scatter(range(len(values)), values, s=2)
            buf = BytesIO()
            plt.savefig(buf, format='jpg', dpi=300)
            buf.seek(0)
            img = Image.open(buf)
            imgs.append(img)
            plt.close()
        imgs[0].save(f'{args.prefix}.{gene}.hap.pdf', "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])

    @classmethod
    def analysis(cls, args):
        for gene in args.gene.split(','):
            ht = cls(args.vcf, args.target, gene, args.proband)
            stat_df = ht.parsing(args)
            haplotype_df = ht.haplotyping(args, stat_df)
            ht.scoring(args, haplotype_df)


def main():
    parser = ArgumentParser()
    parser.add_argument('-target', help='target gene file', required=True)
    parser.add_argument('-gene', help='gene names', required=True)
    parser.add_argument('-vcf', help='vcf file', required=True)
    parser.add_argument('-prefix', help='output file prefix', required=True)
    parser.add_argument('-proband', help='which embryo to proband', default=None)
    parser.add_argument('-correct', help='use correct GT or not? both/neither/trios/embryo', required=True)
    parser.add_argument('-keep_indel', help='use indel or not', action='store_true')
    parser.add_argument('-method_trios', help='ratio adjust trios GT method', default='ratio')
    parser.add_argument('-method_embryo', help='ratio adjust embryo GT method', default='simple')
    parser.add_argument('-triosDepth', help='trios Depth threshold', default=10, type=int)
    parser.add_argument('-embryoDepth', help='embryo Depth threshold', default=10, type=int)
    parser.add_argument('-ref_het_ratio', help='variant ref het ratio threshold', default=0.1, type=float)
    parser.add_argument('-het_hom_ratio', help='variant het hom ratio threshold', default=0.9, type=float)
    parsed_args = parser.parse_args()
    handler = logging.FileHandler(f'{parsed_args.prefix}.log')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    Haplotype.analysis(parsed_args)


if __name__ == '__main__':
    main()
