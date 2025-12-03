from argparse import ArgumentParser
from collections import defaultdict
from io import BytesIO, StringIO
from PIL import Image
import pandas as pd
import subprocess
import logging
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)


# 未指定绝对路径,则环境中需要bcftools
bcftools = 'bcftools'


def get_gene_ps(lfr, chromosome, start, stop):
    # 获取基因上的PS
    while True:
        cmd = f"{bcftools} query -r {chromosome}:{start}-{stop} -f '%CHROM\t%POS\t[%PS\t%GT]\n' {lfr}"
        logger.info(cmd)
        status, output = subprocess.getstatusoutput(cmd)
        if status != 0:
            logger.info('bcftools query error!')
            sys.exit(1)
        block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'PS', 'GT'], dtype={'PS': str, 'GT': str})
        block = block.loc[block['GT'].str.contains('\|') & (block['PS'] != '.')]
        ps = block['PS'].unique().tolist()
        if len(ps) > 0:
            break
        else:
            start -= 1000
            stop += 1000
            if start <= 0:
                start = 1
    return ps


def parse_lfr(lfr):
    # 解析lfr变异
    cmd = f"{bcftools} query -f '%CHROM\t%POS\t[%PS\t%GT]\n' {lfr}"
    logger.info(cmd)
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        logger.info('bcftools query error!')
        sys.exit(1)
    block = pd.read_csv(StringIO(output), sep='\t', header=None, names=['CHROM', 'POS', 'PS', 'GT'], dtype={'PS': str, 'GT': str})
    return block.loc[block['GT'].str.contains('\|') & (block['PS'] != '.')]


def ps_plot(dic, region_dic, out):
    imgs = list()
    for gene in dic:
        if dic[gene]['PS']:
            df = dic[gene]['block'].copy()
            df['value'] = 1
            df['index'] = range(1, df.shape[0] + 1)
            fig, ax = plt.subplots(figsize=(18, 5))
            ax.set_title(f'{gene}', {'fontweight': 'bold'})
            start = df.loc[df['POS'] <= region_dic[gene]['start']].shape[0]
            end = df.loc[df['POS'] <= region_dic[gene]['stop']].shape[0]
            ax.axvline(x=start, ls=':', c='black')
            ax.axvline(x=end, ls=':', c='black')
            for ps in df['PS'].unique():
                df_ps = df.loc[df['PS'] == ps]
                ax.scatter(df_ps['index'].values, df_ps['value'].values, label=ps, s=2)
            ax.legend(title='PS', bbox_to_anchor=(1.05, 1.0), loc='upper left')
            ax.set_xlabel('index')
            ax.axes.yaxis.set_visible(False)
            plt.tight_layout()
            buf = BytesIO()
            plt.savefig(buf, format='jpg', dpi=300)
            buf.seek(0)
            img = Image.open(buf)
            imgs.append(img)
            plt.close()
        else:
            logger.info(f'no PS find in {gene}')
            continue
    if len(imgs) > 0:
        imgs[0].save(out, "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])


def main():
    parser = ArgumentParser()
    parser.add_argument('-target', help='gene region file')
    parser.add_argument('-gene', help='gene list file')
    parser.add_argument('-lfr', help='lfr phased_variants.vcf.gz')
    parser.add_argument('-out', help='output file')
    parsed_args = parser.parse_args()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    region = pd.read_csv(parsed_args.target, sep='\t', header=None, names=['gene', 'chromosome', 'start_extend', 'stop_extend', 'start', 'stop'])
    region.index = region['gene'].tolist()
    region_dic = region.to_dict('index')
    block = parse_lfr(parsed_args.lfr)
    dic = defaultdict(dict)
    genes = pd.read_csv(parsed_args.gene, sep='\t', header=None)
    for gene in genes[0].values:
        ps = get_gene_ps(parsed_args.lfr, region_dic[gene]['chromosome'], region_dic[gene]['start'], region_dic[gene]['stop'])
        dic[gene]['PS'] = ps
        dic[gene]['block'] = block.loc[block['PS'].isin(ps)].copy()
    ps_plot(dic, region_dic, parsed_args.out)


if __name__ == '__main__':
    main()
