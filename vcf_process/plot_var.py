import pandas as pd
from io import StringIO
import subprocess
import sys
import matplotlib.pyplot as plt
import seaborn as sns


def read_vcf(vcf_file):
    grep = 'zgrep' if vcf_file.endswith('.gz') else 'grep'
    com = f'{grep} -v ^## {vcf_file}'
    df = pd.read_csv(StringIO(subprocess.getoutput(com)), sep='\t', dtype=str)
    df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
    df['POS'] = df['POS'].astype(int)
    sample = df.columns[-1]
    df['GT'] = df[sample].str.split(':').str[0]
    gt_list = ['0/1', '1/1', '1/0', '0|1', '1|1', '1|0']
    df = df[df['GT'].isin(gt_list)].copy()
    df['Genotype'] = 'Het'
    df.loc[df['GT'].isin(['1|1', '1/1']), 'Genotype'] = 'Hom'
    df['AD'] = df[sample].str.split(':').str[1]
    df['AD0'] = df['AD'].str.split(',').str[0]
    df['AD1'] = df['AD'].str.split(',').str[1]
    df['Depth'] = df['AD0'].astype(int) + df['AD1'].astype(int)
    df['Ratio'] = (df['AD1'].astype(int) / df['Depth'].astype(int)).round(4)
    df['VarType'] = 'snv'
    df.loc[(df['REF'].str.len() > 1) | (df['ALT'].str.len() > 1), 'VarType'] = 'indel'
    columns = ['CHROM', 'POS', 'REF', 'ALT', 'Depth', 'Ratio', 'VarType', 'Genotype']
    df = df[columns].copy()
    return df


def plot_ratio(df, out_file):
    fig, ax = plt.subplots()
    sns.scatterplot(x='POS', y='Ratio', hue='Genotype', data=df, style='VarType', size='Depth', hue_order=['Het', 'Hom'], style_order=['snv', 'indel'], ax=ax)
    lg = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.savefig(out_file, bbox_extra_artists=(lg,), bbox_inches="tight")
    plt.close(fig)


def main():
    df = read_vcf(sys.argv[1])
    prefix = sys.argv[2]
    regions = sys.argv[3].split(',')
    for region in regions:
        chrom, start, stop = region.split('-')
        _df = df[(df['CHROM'] == chrom) & (df['POS'] >= int(start)) & (df['POS'] <= int(stop))]
        plot_ratio(_df, f'{prefix}.{region}.pdf')


if __name__ == '__main__':
    main()

