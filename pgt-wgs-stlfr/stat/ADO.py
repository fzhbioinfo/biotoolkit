# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
import vcf


class ADO:
    def __init__(self, vcf_file):
        self.reader = vcf.Reader(filename=vcf_file)
        self.embryos = [sample for sample in self.reader.samples if sample.startswith('Embryo')]

    def cal_ado(self, ps_bed_dir, gene, dp_threshold, embryo_dp_threshold, use_indel=False):
        # 计算全基因组和基因区域的ADO,默认不用indel
        if gene == 'all_gene':
            reader = self.reader
        else:
            gene_ps = pd.read_csv(f'{ps_bed_dir}/{gene}.PS_merge.bed', sep='\t', header=None, names=['chromosome', 'start', 'stop'])
            chromosome = gene_ps['chromosome'].values[0]
            start = gene_ps['start'].min()
            stop = gene_ps['stop'].max()
            reader = self.reader.fetch(chromosome, start, stop)
        stat = defaultdict(dict)
        stat['Father']['gene'], stat['Father']['meet'], stat['Father']['ADO'] = gene, '.', '.'
        stat['Mother']['gene'], stat['Mother']['meet'], stat['Mother']['ADO'] = gene, '.', '.'
        if self.embryos:
            for embryo in self.embryos:
                stat[embryo]['gene'], stat[embryo]['meet'], stat[embryo]['ADO'] = gene, 0, 0
            for record in reader:
                if gene == 'all_gene' and (record.CHROM in ['chrX', 'chrY'] or 'M' in record.CHROM):
                    continue
                try:
                    if 'LowQual' in record.FILTER:
                        continue
                except TypeError:
                    pass
                if '*' in record.ALT:
                    continue
                if not use_indel:
                    if record.var_type != 'snp':
                        continue
                if record.genotype('Father')['GT'] in ['./.', '.|.'] or record.genotype('Mother')['GT'] in ['./.', '.|.']:
                    continue
                try:
                    ad_dad = record.genotype('Father')['AD']
                    ad_mom = record.genotype('Mother')['AD']
                    dp_dad = sum(ad_dad)
                    dp_mom = sum(ad_mom)
                    gt_dad = record.genotype('Father')['GT'].replace('|', '/')
                    gt_mom = record.genotype('Mother')['GT'].replace('|', '/')
                    if dp_dad < dp_threshold or dp_mom < dp_threshold:
                        continue
                    allele_dad = set(gt_dad.split('/'))
                    allele_mom = set(gt_mom.split('/'))
                    if len(allele_dad) != 1 or len(allele_mom) != 1 or allele_dad == allele_mom:
                        continue
                except AttributeError:
                    continue
                for embryo in self.embryos:
                    stat[embryo]['meet'] += 1
                    gt = record.genotype(embryo)['GT'].replace('|', '/')
                    try:
                        ad = record.genotype(embryo)['AD']
                        dp = sum(ad)
                    except AttributeError:
                        continue
                    if len(set(gt.split('/'))) == 1 and gt.split('/')[0] in [gt_dad[0], gt_mom[0]] and dp >= embryo_dp_threshold:
                        stat[embryo]['ADO'] += 1
        stat_df = pd.DataFrame.from_dict(stat, orient='index')
        return stat_df

    @classmethod
    def analysis(cls, args):
        genes = pd.read_csv(args.gene, sep='\t', header=None, names=['gene'])
        ado = cls(args.family_vcf)
        ado_stat = pd.DataFrame()
        for gene in ['all_gene'] + genes['gene'].tolist():
            stat_df = ado.cal_ado(args.ps_bed_dir, gene, args.triosDepth, args.embryoDepth)
            ado_stat = ado_stat.append(stat_df, sort=False)
        ado_stat.to_csv(args.out, sep='\t')


def main():
    parser = ArgumentParser()
    parser.add_argument('-family_vcf', help='vcf file')
    parser.add_argument('-ps_bed_dir', help='ps bed dir')
    parser.add_argument('-gene', help='gene file')
    parser.add_argument('-use_indel', help='use indel or not', action='store_true')
    parser.add_argument('-triosDepth', help='depth threshold', default=10, type=int)
    parser.add_argument('-embryoDepth', help='embryo depth threshold', default=4, type=int)
    parser.add_argument('-out', help='out file', required=True)
    parsed_args = parser.parse_args()
    ADO.analysis(parsed_args)


if __name__ == '__main__':
    main()
