from argparse import ArgumentParser
from io import BytesIO
import pandas as pd
import numpy as np
import os
import allel
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
mpl.rcParams['hatch.color'] = 'black'
ROOT = os.path.abspath(os.path.dirname(__file__))


class Ideogram:
    def __init__(self, cytoband, gap):
        self.cytoband = pd.read_csv(cytoband, sep='\t', header=None, names=['chromosome', 'start', 'stop', 'band', 'stain'])
        self.cytoband['gray'] = self.cytoband['stain'].str.extract(r'(\d+)$')
        self.cytoband['gray'] = self.cytoband['gray'].fillna(0).astype(int) / 100
        self.gap = pd.read_csv(gap, sep='\t')

    def plot(self, chromosome, ymin, ymax, ratio, column):
        cytoband = self.cytoband.loc[self.cytoband['chromosome'] == chromosome].copy()
        end = cytoband['stop'].max()
        cytoband['start'] = cytoband['start'] / end
        cytoband['stop'] = cytoband['stop'] / end
        gap = self.gap.loc[self.gap['chromosome'] == chromosome].copy()
        gap['start'] = gap['start'] / end
        gap['stop'] = gap['stop'] / end
        fig, ax = plt.subplots(2, 1, figsize=(18, 6), gridspec_kw={'height_ratios': [0.95, 0.05]})
        # cnv散点图
        ax[0].set_ylabel('Ratio')
        het = ratio.loc[ratio['genotype'] == 'het']
        hom = ratio.loc[ratio['genotype'] == 'hom']
        # ax[0].scatter(ratio['pos'].values / end, ratio[column].values, c='black', s=0.1)
        ax[0].scatter(het['pos'].values / end, het[column].values, c='blueviolet', s=0.1)
        ax[0].scatter(hom['pos'].values / end, hom[column].values, c='orange', s=0.1)
        ax[0].set_xlim(0, 1)
        ax[0].set_ylim(ymin, ymax)
        ax[0].set_xticks([])
        ax[0].set_title(chromosome)
        for row in gap.itertuples():
            if row.gap_type == 'pseudo_autosomal':
                ax[0].fill_between([row.start, row.stop], ymin, ymax, color='green', alpha=0.8)
            else:
                ax[0].fill_between([row.start, row.stop], ymin, ymax, color='gray', alpha=0.8)
        # 染色体G带图
        i = 1
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(-1.2, 1.2)
        for row in cytoband.itertuples():
            if row.stain == 'acen':
                x = np.arange(row.start, row.stop, (row.stop - row.start) / 1000)
                if i == 1:
                    k = 1 / (row.start - row.stop)
                    y1 = k * (x - row.stop)
                    y2 = k * (row.stop - x)
                    ax[1].fill_between(x, y2, y1, color='red')
                    ax[1].axhline(y=1, xmin=0, xmax=row.start, color='gray')
                    ax[1].axhline(y=-1, xmin=0, xmax=row.start, color='gray')
                    i += 1
                else:
                    k = 1 / (row.stop - row.start)
                    y1 = k * (x - row.start)
                    y2 = k * (row.start - x)
                    ax[1].fill_between(x, y1, y2, color='red')
                    ax[1].axhline(y=1, xmin=row.stop, xmax=1, color='gray')
                    ax[1].axhline(y=-1, xmin=row.stop, xmax=1, color='gray')
            elif row.stain == 'gvar':
                ax[1].axhspan(-1., 1., xmin=row.start, xmax=row.stop, hatch='//')
            else:
                ax[1].axhspan(-1., 1., xmin=row.start, xmax=row.stop, color='black', alpha=row.gray)
            ax[1].text((row.start + row.stop) / 2, -1.1, row.band, rotation='vertical', verticalalignment='top')
        ax[1].axvline(x=0, ymin=-1, ymax=1, color='gray')
        ax[1].axvline(x=1, ymin=-1, ymax=1, color='gray')
        plt.axis('off')
        plt.tight_layout()

    @classmethod
    def run(cls, args):
        imgs = list()
        ide = cls(args.cyto, args.gap)
        # ratio格式处理
        ratio = parse_vcf(args.vcf)
        ratio = ratio.loc[~ratio['chromosome'].str.contains('chrM')].copy()
        ratio['order'] = ratio['chromosome'].str.extract(r'chr(.*?)$')
        ratio.loc[ratio['order'] == 'X', 'order'] = '23'
        ratio.loc[ratio['order'] == 'Y', 'order'] = '24'
        ratio['order'] = ratio['order'].astype(int)
        ratio.sort_values(by=['order'], inplace=True)
        if 'pos' not in ratio.columns:
            ratio['pos'] = (ratio['start'] + ratio['stop']) / 2
        # 按染色体画图
        chromosomes = [c for c in args.chromosome.split(',')] if args.chromosome else ratio['chromosome'].unique()
        for chromosome in chromosomes:
            df = ratio.loc[ratio['chromosome'] == chromosome]
            ide.plot(chromosome, args.ymin, args.ymax, df, args.col)
            buf = BytesIO()
            plt.savefig(buf, format='jpg', dpi=300)
            buf.seek(0)
            img = Image.open(buf)
            imgs.append(img)
            plt.close()
        imgs[0].save(args.out, "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])


def parse_vcf(vcf):
    df = pd.DataFrame()
    dataset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 'calldata/AD', 'variants/is_snp',
                                          'calldata/DP', 'calldata/GT'])
    depth_cut = dataset['calldata/DP'].flatten() >= 10
    is_snp = dataset['variants/is_snp']
    df['chromosome'] = dataset['variants/CHROM'][depth_cut & is_snp]
    df['pos'] = dataset['variants/POS'][depth_cut & is_snp]
    df['DP'] = dataset['calldata/DP'].squeeze()[depth_cut & is_snp]
    df['DP_alt'] = dataset['calldata/AD'].squeeze()[:, 1][depth_cut & is_snp]
    df['ratio'] = df['DP_alt'] / df['DP']
    df['genotype'] = 'hom'
    df.loc[allel.GenotypeArray(dataset['calldata/GT'])[depth_cut & is_snp].is_het().flatten(), 'genotype'] = 'het'
    return df


def main():
    parser = ArgumentParser()
    parser.add_argument('-cyto', help='cytoBand', default=f'{ROOT}/cytoBand.txt.gz')
    parser.add_argument('-gap', help='chromosome gap region', default=f'{ROOT}/GRCh37.gap.regions.tsv')
    parser.add_argument('-ymin', help='y min', default=0, type=float)
    parser.add_argument('-ymax', help='y max', default=4, type=float)
    parser.add_argument('-ratio', help='copy ratio file')
    parser.add_argument('-out', help='output file', required=True)
    parser.add_argument('-col', help='column to plot', default='CopyRatio')
    parser.add_argument('-chromosome', help='which chromosome to plot')
    parser.add_argument('-vcf', help='vcf file')
    parsed_args = parser.parse_args()
    Ideogram.run(parsed_args)


if __name__ == '__main__':
    main()
