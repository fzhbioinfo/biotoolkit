# -*- coding:utf-8 -*-
from collections import defaultdict
from argparse import ArgumentParser
import pandas as pd
import vcf
import re


class Annotation:
    def __init__(self, vcf_file):
        self.reader = vcf.Reader(filename=vcf_file)

    def parsing(self):
        try:
            ann = [i.replace(' ', '') for i in re.findall('\'(.*?)\'', self.reader.infos['ANN'].desc)[0].split('|')]
        except KeyError:
            pass
        annotation = list()
        for record in self.reader:
            if '*' in record.ALT:
                continue
            for info in record.INFO['ANN']:
                allele = info.split('|')[ann.index('Allele')]
                func = info.split('|')[ann.index('Annotation')]
                gene = info.split('|')[ann.index('Gene_Name')]
                transcript = info.split('|')[ann.index('Feature_ID')]
                hgvs_c = info.split('|')[ann.index('HGVS.c')]
                hgvs_p = info.split('|')[ann.index('HGVS.p')]
                for sample in self.reader.samples:
                    call = defaultdict(dict)
                    call['CHROM'] = record.CHROM
                    call['POS'] = record.POS
                    call['REF'] = record.REF
                    alts = [str(alt) for alt in record.ALT]
                    call['ALT'] = ','.join(alts)
                    call['VarType'] = record.var_type
                    call['sample'] = sample
                    gt = record.genotype(sample)['GT'].replace('|', '/')
                    call['GT'] = gt
                    try:
                        if gt == './.':
                            raise ValueError
                        ad = record.genotype(sample)['AD']
                        dp = sum(ad)
                        call['DP'] = dp
                        call['AD'] = ','.join([str(i) for i in ad])
                    except (AttributeError, ValueError):
                        call['DP'] = '.'
                        call['AD'] = '.'
                    call['Allele'] = allele
                    call['Function'] = func
                    call['Gene Symbol'] = gene
                    call['Transcript'] = transcript
                    call['HGVS.c'] = hgvs_c
                    call['HGVS.p'] = hgvs_p
                    annotation.append(call)
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'VarType', 'sample', 'GT', 'DP', 'AD', 'Allele', 'Function', 'Gene Symbol', 'Transcript', 'HGVS.c', 'HGVS.p']
        annotation_df = pd.DataFrame(annotation, columns=columns)
        return annotation_df

    @classmethod
    def analysis(cls, args):
        at = cls(args.vcf)
        df = at.parsing()
        df.to_csv(args.out, sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-vcf', help='vcf file', required=True)
    parser.add_argument('-out', help='output file', required=True)
    parsed_args = parser.parse_args()
    Annotation.analysis(parsed_args)


if __name__ == '__main__':
    main()
