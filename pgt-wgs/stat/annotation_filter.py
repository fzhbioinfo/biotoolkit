# -*- coding:utf-8 -*-
from argparse import ArgumentParser
import pandas as pd
import os


def annotation_filter(args):
    genes = [g for g in args.gene.split(',')]
    gene_trans = pd.read_csv(args.trans, sep='\t', usecols=['gene', 'Gene Symbol', 'Transcript'])
    annotation = pd.read_csv(args.anno, sep='\t')
    dipin = pd.read_csv(args.dipin, sep='\t')
    annotation = pd.merge(annotation, dipin[['Common Name', 'Transcript', 'HGVS.c']], on=['Transcript', 'HGVS.c'], how='left')
    trio = os.path.basename(os.path.dirname(os.path.dirname(args.anno)))
    info = pd.read_csv(args.info, sep='\t', usecols=['sample', 'sampleID'])
    info.drop_duplicates(subset=['sampleID'], inplace=True)
    info['ID'] = info['sample'] + '_' + info['sampleID']
    info = info.loc[(info['sample'].isin(['Father', 'Mother'])) | (info['ID'] == trio) | (info['sample'].str.contains('Embryo'))].copy()
    info.drop(columns=['ID'], inplace=True)
    annotation = pd.merge(annotation, info, on=['sample'])
    annotation.insert(6, 'sampleID', annotation.pop('sampleID'))
    annotation.sort_values(by=['CHROM', 'POS'], inplace=True)
    annotation.fillna('.', inplace=True)
    for gene in genes:
        if gene != 'HLA':
            gene_trans_ = gene_trans[gene_trans['gene'] == gene]
            df = annotation[(annotation['Gene Symbol'].isin(gene_trans_['Gene Symbol'].tolist())) & (annotation['Transcript'].isin(gene_trans_['Transcript'].tolist()))]
            df.to_csv(f'{args.prefix}.{gene}.anno.tsv', sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-gene', help='genes', required=True)
    parser.add_argument('-trans', help='gene trans', required=True)
    parser.add_argument('-anno', help='annotation result', required=True)
    parser.add_argument('-prefix', help='output file prefix', required=True)
    parser.add_argument('-info', help='sample info', required=True)
    parser.add_argument('-dipin', help='dipin common name file', required=True)
    parsed_args = parser.parse_args()
    annotation_filter(parsed_args)


if __name__ == '__main__':
    main()
