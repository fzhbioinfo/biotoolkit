# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
from io import BytesIO, StringIO
from interval import Interval
from PIL import Image
import pandas as pd
import numpy as np
import subprocess
import logging
import sys
import vcf
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)

# 未指定绝对路径,则环境中需要bcftools
bcftools = 'bcftools'


class Haplotype:
    def __init__(self, family_vcf):
        # 类初始化,读取家系vcf,需包括Mother和Father,获取胚胎列表
        self.family_vcf = family_vcf
        self.reader = vcf.Reader(filename=self.family_vcf)
        self.embryos = [int(sample.split('Embryo')[1]) for sample in self.reader.samples if sample.startswith('Embryo')]
        self.embryos.sort()
        self.embryos = [f'Embryo{i}' for i in self.embryos]
        self.family = ['Mother', 'Father']

    @staticmethod
    def parse_cnv_coordinate(cnv):
        # 解析cnv染色体坐标
        chromosome, start, stop = re.findall('(.*?)[:-](.*?)[:-](.*?)$', cnv)[0]
        chromosome = f'chr{chromosome}' if not chromosome.startswith('chr') else chromosome
        return chromosome, int(start), int(stop)

    @staticmethod
    def parse_snv_coordinate(snv):
        # 解析snp/indel染色体坐标
        chromosome, pos, ref, alt = re.findall('(.*?)[:-](.*?)[:-](.*?)[:-](.*?)$', snv)[0]
        chromosome = f'chr{chromosome}' if not chromosome.startswith('chr') else chromosome
        return chromosome, int(pos), ref, alt

    @staticmethod
    def get_cnv_ps(lfr, chromosome, start, stop):
        # 从LFR的cnv分型结果检索目标cnv,若存在多个则取交集最大的cnv,并获取其PS和判断在哪条allele
        cnv = Interval(start, stop)
        if stop - start <= 30000:
            reader = vcf.Reader(filename=f'{lfr}/outs/dels.vcf.gz')
        else:
            reader = vcf.Reader(filename=f'{lfr}/outs/large_svs.vcf.gz')
        sample = reader.samples[0]
        records = reader.fetch(chromosome, start, stop)
        # 默认为None
        ps, allele_number = None, None
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
                allele_number = judge_cnv_ps_allele(record.genotype(sample)['GT'])
        return ps, allele_number

    @staticmethod
    def get_snv_ps(lfr, chromosome, pos, ref, alt):
        # 读取LFR的snp/indel变异分型结果,提取目标变异的PS和判断在哪条allele
        reader = vcf.Reader(filename=f'{lfr}/outs/phased_variants.vcf.gz')
        sample = reader.samples[0]
        records = reader.fetch(chromosome, pos - 1, pos)
        # 默认为None
        ps, allele_number = None, None
        for record in records:
            try:
                if record.genotype(sample)['PS'] is None or '/' in record.genotype(sample)['GT']:
                    continue
            except KeyError:
                continue
            alt_lfr = [str(a) for a in record.ALT]
            if chromosome == record.CHROM and pos == record.POS and ref == record.REF and alt in alt_lfr:
                ps = record.genotype(sample)['PS']
                allele_number = judge_snv_ps_allele(record.genotype(sample)['GT'], alt_lfr, alt)
                break
        return ps, allele_number

    @staticmethod
    def parse_ps_block(lfr, chromosome, ps):
        # bcftools解析LFR的变异分型结果,获取目标PS的block和坐标范围,并过滤无效变异
        cmd = f"{bcftools} query -r {chromosome} -f '%CHROM\t%POS\t%REF\t[%PS\t%GT\t%TGT\t%DP\t%AD]\n' {lfr}/outs/phased_variants.vcf.gz"
        logger.info(cmd)
        status, output = subprocess.getstatusoutput(cmd)
        if status != 0:
            logger.info('bcftools query error!')
            sys.exit(1)
        block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'PS', 'GT_lfr', 'TGT_lfr', 'DP_lfr', 'AD_lfr'],
                            dtype={'PS': str, 'DP_lfr': str})
        block = block[(block['PS'] == str(ps)) & (block.apply(lambda x: x['DP_lfr'] != '.' and int(x['DP_lfr']) > 0, axis=1))].copy()
        block['DP_lfr'] = block['DP_lfr'].astype(int)
        block['Ratio_lfr'] = (block['AD_lfr'].str.split(',').str[-1].astype(int) / block['DP_lfr']).round(3)
        return block, block['POS'].min(), block['POS'].max()

    def parse_family_block(self, chromosome, start, stop):
        # 解析block范围内家系vcf的变异信息
        reader = self.reader.fetch(chromosome, start - 1, stop)
        stat = list()
        for record in reader:
            # 父/母基因型分析失败该变异记录舍弃
            skip = False
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
            for member in self.family:
                if record.genotype(member)['GT'] in ['./.', '.|.']:
                    call = defaultdict(dict)
                    skip = True
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
                    dp = sum(ad)
                    ratio = round(ad[-1] / dp, 3)
                    call[f'{embryo}_GT'] = gt
                    call[f'{embryo}_TGT'] = '/'.join([base[int(i)] for i in gt.split('/')])
                    call[f'{embryo}_DP'] = dp
                    call[f'{embryo}_AD'] = ','.join([str(i) for i in ad])
                    call[f'{embryo}_Ratio'] = ratio
                except (AttributeError, ZeroDivisionError, ValueError):
                    call[f'{embryo}_GT'] = '.'
                    call[f'{embryo}_TGT'] = '.'
                    call[f'{embryo}_DP'] = '.'
                    call[f'{embryo}_AD'] = '.'
                    call[f'{embryo}_Ratio'] = '.'
            if call:
                stat.append(call)
        stat_df = pd.DataFrame(stat)
        return stat_df

    def haplotyping(self, family_block_lfr, tag, args):
        # LFR分型来源于Maternal,则另一方为Father;LFR分型来源于Paternal,则另一方为Mother
        other = 'Father' if tag == 'M' else 'Mother'
        result = list()
        # 迭代每条变异记录
        for record in family_block_lfr.itertuples(index=False):
            # 转为OrderDict处理简便
            dic = record._asdict()
            # 父母的两条单体碱基序列
            tgt_lfr = dic['TGT_lfr'].split('|')
            tgt_other = dic[f'{other}_TGT'].split('/')
            # IF判断,LFR分型结果需为杂合,父母相同的碱基序列小于2
            if len(set(tgt_lfr)) <= 1 or dic[f'{other}_TGT'] in ['.', './.']:
                dic['IF'] = 'no'
            elif len(set(tgt_other)) == 1 or len(set(tgt_other) & set(tgt_lfr)) < 2:
                dic['IF'] = 'yes'
            else:
                dic['IF'] = 'no'
            # IF可用来判断哪条allele,是否通过深度过滤
            if dic['IF'] == 'no':
                dic['Allele'] = '.'
                dic['IF_filter'] = 'Fail'
            else:
                dic['Allele'] = self.which_allele(tgt_other, [dic[f'{tag}1'], dic[f'{tag}2']], tag)
                dic['IF_filter'] = 'LowDepth' if dic[f'{other}_DP'] < args.triosDepth or dic['DP_lfr'] < args.triosDepth else 'Pass'
            # 根据胚胎的碱基序列判断其遗传的单体,深度低的变异加上LowDepth标签
            for embryo in self.embryos:
                if dic['IF'] == 'no' or dic[f'{embryo}_TGT'] == '.':
                    dic[f'{embryo}_hap'] = '.'
                else:
                    # 是IF的位点才进行遗传单体判断
                    dic[f'{embryo}_hap'] = self.support_allele(tgt_other, [dic[f'{tag}1'], dic[f'{tag}2']], dic[f'{embryo}_TGT'].split('/'), tag)
                    if dic[f'{embryo}_DP'] < args.embryoDepth and dic[f'{embryo}_hap'] != '.':
                        dic[f'{embryo}_hap'] = dic[f'{embryo}_hap'] + '-LowDepth'
            result.append(dic)
        haplotype = pd.DataFrame(result)
        return haplotype

    @staticmethod
    def which_allele(tgt_other, hap_lfr, tag):
        # IF可用于区分父/母哪条单体,父母共同的碱基数<=1
        common = list(set(tgt_other) & set(hap_lfr))
        # 可以区分LFR特有碱基所在的单体
        allele = f'{tag}{2 - hap_lfr.index(common[0])}' if common else f'{tag}1|{tag}2'
        return allele

    @staticmethod
    def support_allele(tgt_other, hap_lfr, tgt_embryo, tag):
        # 胚胎纯合
        if len(set(tgt_embryo)) == 1:
            # 胚胎的碱基序列存在于亲代other一方或与亲代LFR一方无交集则无法进行遗传单体判断,反之则遗传与亲代LFR一方共同碱基所在的单体
            if set(tgt_embryo) & set(tgt_other) or not set(tgt_embryo) & set(hap_lfr):
                hap = '.'
            else:
                hap = f'{tag}{hap_lfr.index(tgt_embryo[0]) + 1}'
        # 胚胎杂合
        else:
            # 胚胎与亲代LFR一方相同的碱基
            common = list(set(tgt_embryo) & set(hap_lfr))
            # 一个相同碱基且不在亲代other一方,则遗传亲代LFR一方此碱基所在的单体
            if len(common) == 1 and common[0] not in tgt_other:
                hap = f'{tag}{hap_lfr.index(common[0]) + 1}'
            # 两个相同碱基且亲代other一方为纯合,同时父母有相同碱基序列,则遗传亲代LFR一方特有碱基所在的单体
            elif len(common) == 2 and len(set(tgt_other)) == 1 and tgt_other[0] in hap_lfr:
                common.remove(tgt_other[0])
                hap = f'{tag}{hap_lfr.index(common[0]) + 1}'
            else:
                hap = '.'
        return hap

    @staticmethod
    def counting(upstream, internal, downstream, sample, variant):
        # 支持的单体型计数
        count = defaultdict(dict)
        count['sample'] = sample
        count['Variant'] = variant
        column = 'Allele' if sample == 'IF' else f'{sample}_hap'
        upstream_count = upstream[column].value_counts().to_dict()
        internal_count = internal[column].value_counts().to_dict()
        downstream_count = downstream[column].value_counts().to_dict()
        for i in ['F1', 'F2', 'M1', 'M2']:
            if sample == 'IF':
                hap12 = 'F1|F2' if i in ['F1', 'F2'] else 'M1|M2'
                upstream_count[i] = upstream_count.get(i, 0) + upstream_count.get(hap12, 0)
                internal_count[i] = internal_count.get(i, 0) + internal_count.get(hap12, 0)
                downstream_count[i] = downstream_count.get(i, 0) + downstream_count.get(hap12, 0)
            count[f'{i}-U'] = upstream_count.get(i, 0)
            count[f'{i}-I'] = internal_count.get(i, 0)
            count[f'{i}-D'] = downstream_count.get(i, 0)
        return count

    @staticmethod
    def resulting(score_df, embryos, up=0.2, low=0.05):
        # 单体型判断结果,目标单体型IF比例0.2,非目标0.05
        score_df.index = score_df['sample'].tolist()
        for embryo in embryos:
            for tag in ['F', 'M']:
                result = list()
                for uid in ['U', 'I', 'D']:
                    if_1st = score_df.loc['IF', f'{tag}1-{uid}']
                    if_2nd = score_df.loc['IF', f'{tag}2-{uid}']
                    embryo_1st = score_df.loc[embryo, f'{tag}1-{uid}']
                    embryo_2nd = score_df.loc[embryo, f'{tag}2-{uid}']
                    if if_1st > 0 and if_2nd > 0:
                        r1 = embryo_1st / if_1st
                        r2 = embryo_2nd / if_2nd
                        if r1 > up and r2 < low:
                            result.append(f'{tag}1')
                        elif r1 < low and r2 > up:
                            result.append(f'{tag}2')
                        else:
                            result.append('-')
                    elif if_1st > 0 and if_2nd == 0:
                        r1 = embryo_1st / if_1st
                        if r1 > up:
                            result.append(f'{tag}1')
                        else:
                            result.append('-')
                    elif if_1st == 0 and if_2nd > 0:
                        r2 = embryo_2nd / if_2nd
                        if r2 > up:
                            result.append(f'{tag}2')
                        else:
                            result.append('-')
                    else:
                        result.append('-')
                if tag == 'F':
                    score_df.loc[embryo, 'Paternal'] = '|'.join(result)
                else:
                    score_df.loc[embryo, 'Maternal'] = '|'.join(result)
        return score_df

    def scoring(self, haplotype_df, parent, args):
        # 结果过滤
        if not args.keep_indel:
            haplotype_df = haplotype_df[haplotype_df['VarType'] != 'indel']
        haplotype_df = haplotype_df[haplotype_df['IF_filter'] == 'Pass']
        # 如果有胚胎样本则画图
        if self.embryos and not haplotype_df.empty:
            self.scatter(haplotype_df, parent, args)
        # 变异维度统计
        score_all = pd.DataFrame(columns=['sample', 'Variant'] + [f'{i}{k}-{j}' for i in 'FM' for k in '12' for j in 'UID' ] + ['Paternal', 'Maternal'])
        for var in haplotype_df['Variant'].unique().tolist():
            score = list()
            var_df = haplotype_df.loc[haplotype_df['Variant'] == var]
            upstream = var_df[var_df['Location'] == 'upstream']
            internal = var_df[var_df['Location'] == 'internal']
            downstream = var_df[var_df['Location'] == 'downstream']
            # IF与胚胎分型计数
            for sample in ['IF'] + self.embryos:
                score.append(self.counting(upstream, internal, downstream, sample, var))
            score_df = pd.DataFrame(score)
            score_df['Paternal'], score_df['Maternal'] = '', ''
            # 存在胚胎样本则进行遗传单体自动判断
            if self.embryos:
                score_df = self.resulting(score_df, self.embryos)
            # 多个变异结果合并
            score_all = score_all.append(score_df)
        score_all.to_excel(f'{args.prefix}.{parent}.score.xlsx', index=False)

    def scatter(self, haplotype_df, parent, args):
        # 单体型结果画图
        imgs = list()
        ticks = ['M1', 'M2'] if parent.startswith('M') else ['F1', 'F2']
        y_locator = plt.MultipleLocator(1)
        for var in haplotype_df['Variant'].unique().tolist():
            var_df = haplotype_df.loc[haplotype_df['Variant'] == var]
            for embryo in self.embryos:
                fig, ax = plt.subplots(figsize=(10, 3))
                ax.yaxis.set_major_locator(y_locator)
                ax.set_ylim(0.5, 2.5)
                plt.yticks(range(1, 3), ticks)
                ax.set_title(f'{embryo}_{var}', {'fontweight': 'bold'})
                hap = var_df[var_df[f'{embryo}_hap'].isin(ticks)]
                start = hap.loc[hap['Location'] == 'upstream'].shape[0]
                end = hap.loc[hap['Location'].isin(['upstream', 'internal'])].shape[0]
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
        imgs[0].save(f'{args.prefix}.{parent}.hap.pdf', "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])

    def block_analysis(self, genes, variant, var_type, lfr, parent, args):
        # 读取基因的坐标信息
        region = pd.read_csv(args.target, sep='\t', header=None, names=['gene', 'chromosome', 'start_extend', 'stop_extend', 'start', 'stop'])
        region.index = region['gene'].tolist()
        region_dic = region.to_dict('index')
        # 判断LFR单体分型来源父/母
        tag = 'M' if parent.startswith('M') else 'F'
        # 目标基因,变异,变异类型
        gene_list = [g for g in genes.split(',')]
        var_list = [v for v in variant.split(',')]
        var_type_list = [t for t in var_type.split(',')]
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'VarType'] + [f'{s}_{i}' for s in self.family + self.embryos for i in ['GT', 'TGT', 'DP', 'AD', 'Ratio']] + \
                  ['PS'] + [f'{i}_lfr' for i in ['GT', 'TGT', 'DP', 'AD', 'Ratio']] + [f'{tag}1', f'{tag}2', 'IF', 'Allele', 'IF_filter'] + [f'{s}_hap' for s in self.embryos] + ['Target', 'Location', 'Variant']
        haplotype_df = pd.DataFrame(columns=columns)
        for i, gene in enumerate(gene_list):
            var = var_list[i]
            if var_type_list[i] == 'cnv':
                chromosome, start, stop = self.parse_cnv_coordinate(var)
                ps, allele_number = self.get_cnv_ps(lfr, chromosome, start, stop)
            elif var_type_list[i] in ['snv', 'indel']:
                chromosome, pos, ref, alt = self.parse_snv_coordinate(var)
                ps, allele_number = self.get_snv_ps(lfr, chromosome, pos, ref, alt)
            else:
                # 待开发...
                chromosome = None
                ps, allele_number = None, None
            if ps is None:
                logger.info(f'Can not find effective PS! Haplotyping {parent} {var_type_list[i]}:{var} failed!')
            else:
                block, pos_min, pos_max = self.parse_ps_block(lfr, chromosome, ps)
                family_block = self.parse_family_block(chromosome, pos_min, pos_max)
                # 合并LFR分型信息
                family_block_lfr = pd.merge(family_block, block, on=['CHROM', 'POS', 'REF'], how='left')
                family_block_lfr = family_block_lfr[family_block_lfr['TGT_lfr'].notna()].copy()
                # 确定父/母两条单体型
                family_block_lfr[f'{tag}1'] = family_block_lfr['TGT_lfr'].str.replace('/', '|').str.split('|').str[allele_number]
                family_block_lfr[f'{tag}2'] = family_block_lfr['TGT_lfr'].str.replace('/', '|').str.split('|').str[1 - allele_number]
                # 胚胎遗传单体分析
                haplotype = self.haplotyping(family_block_lfr, tag, args)
                # 是否是目标变异
                haplotype['Target'] = '.'
                if var_type_list[i] in ['snv', 'indel']:
                    haplotype.loc[(haplotype['CHROM'] == chromosome) & (haplotype['POS'] == pos) & (haplotype['REF'] == ref), 'Target'] = 'yes'
                # 添加每条变异记录与基因的位置关系
                haplotype['Location'] = 'internal'
                haplotype.loc[haplotype['POS'] < region_dic[gene]['start'], 'Location'] = 'upstream'
                haplotype.loc[haplotype['POS'] > region_dic[gene]['stop'], 'Location'] = 'downstream'
                # 添加目标变异名,多个目标分型结果合并
                haplotype['Variant'] = f'{gene}_{var_type_list[i]}_{var}'
                # 对分型结果的变异进行过滤
                if args.extend:
                    haplotype = haplotype.loc[(haplotype['POS'] > region_dic[gene]['start'] - args.left) & (haplotype['POS'] < region_dic[gene]['stop'] + args.right)]
                haplotype_df = haplotype_df.append(haplotype)
        haplotype_df.to_csv(f'{args.prefix}.{parent}.hap.tsv', sep='\t', index=False)
        # 输出统计结果
        self.scoring(haplotype_df, parent, args)

    @classmethod
    def analysis(cls, args):
        ht = cls(args.family_vcf)
        if args.info:
            info = pd.read_csv(args.info, sep='\t')
            info = info.loc[info['is_lfr'] == 'yes'].copy()
            info.drop_duplicates(subset=['sample'], inplace=True)
            info.index = info['sample'].values
            info_dic = info.to_dict('index')
            if 'Mother' in info_dic.keys():
                ht.block_analysis(info_dic['Mother']['gene'], info_dic['Mother']['variant'], info_dic['Mother']['var_type'], f'{args.lfr_dir}/Mother', 'Maternal', args)
            if 'Father' in info_dic.keys():
                ht.block_analysis(info_dic['Father']['gene'], info_dic['Father']['variant'], info_dic['Father']['var_type'], f'{args.lfr_dir}/Father', 'Paternal', args)
        else:
            if args.mat_gene:
                ht.block_analysis(args.mat_gene, args.mat_var, args.mat_var_type, args.mat_lfr, 'Maternal', args)
            if args.pat_gene:
                ht.block_analysis(args.pat_gene, args.pat_var, args.pat_var_type, args.pat_lfr, 'Paternal', args)


def judge_cnv_ps_allele(gt):
    # 根据LFR分型结果判断cnv在哪条allele
    if gt in ['1|0', '1|1']:
        allele_number = 0
    elif gt == '0|1':
        allele_number = 1
    else:
        allele_number = None
    return allele_number


def judge_snv_ps_allele(gt, alt_lfr, alt):
    # 根据LFR分型结果判断snv在哪条allele
    if gt in ['1|0', '1|1']:
        allele_number = 0
    elif gt == '0|1':
        allele_number = 1
    elif gt == '1|2':
        allele_number = alt_lfr.index(alt)
    elif gt == '2|1':
        allele_number = 1 - alt_lfr.index(alt)
    else:
        allele_number = None
    return allele_number


def main():
    parser = ArgumentParser()
    parser.add_argument('-target', help='gene region file')
    parser.add_argument('-pat_gene', help='paternal genes', default=None)
    parser.add_argument('-mat_gene', help='maternal genes', default=None)
    parser.add_argument('-pat_var', help='paternal variant', default=None)
    parser.add_argument('-mat_var', help='maternal variant', default=None)
    parser.add_argument('-pat_var_type', help='paternal variant type', default=None)
    parser.add_argument('-mat_var_type', help='maternal variant type', default=None)
    parser.add_argument('-pat_lfr', help='paternal stLFR dir', default=None)
    parser.add_argument('-mat_lfr', help='maternal stLFR dir', default=None)
    parser.add_argument('-family_vcf', help='family vcf file', required=True)
    parser.add_argument('-prefix', help='output file prefix', required=True)
    parser.add_argument('-keep_indel', help='use indel or not', action='store_true')
    parser.add_argument('-triosDepth', help='trios Depth threshold', default=10, type=int)
    parser.add_argument('-embryoDepth', help='embryo Depth threshold', default=10, type=int)
    parser.add_argument('-left', help='gene left extend', default=2000000, type=int)
    parser.add_argument('-right', help='gene right extend', default=2000000, type=int)
    parser.add_argument('-extend', help='gene extend or not', action='store_true')
    parser.add_argument('-info', help='sample info', default=None)
    parser.add_argument('-lfr_dir', help='lfr dir', default=None)
    parsed_args = parser.parse_args()
    handler = logging.FileHandler(f'{parsed_args.prefix}.log')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    Haplotype.analysis(parsed_args)


if __name__ == '__main__':
    main()
