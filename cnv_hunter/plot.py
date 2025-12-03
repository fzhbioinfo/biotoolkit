"""画图
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_segment(df, cnv_name, column, out_file):
    """画出CNV片段
    """
    sample = df['sample'].values[0]
    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    for i, chrom in enumerate(['chr11', 'chr16']):
        df_chrom = df[df['chrom'] == chrom].copy()
        df_chrom.index = range(df_chrom.shape[0])
        df_chrom_normal = df_chrom[df_chrom['type'] == 'NORMAL']
        df_chrom_cnv = df_chrom[df_chrom['type'] != 'NORMAL']
        ax[i].set_ylabel('copy ratio')
        ax[i].set_xlabel('bin number')
        ax[i].set_title(f'{sample} {chrom}')
        ax[i].set_ylim(0, 5)
        ax[i].scatter(df_chrom_normal.index, df_chrom_normal[column] * 2, color='green')
        ax[i].scatter(df_chrom_cnv.index, df_chrom_cnv[column] * 2, color='red')
        y = 1
        for row in cnv_name.loc[cnv_name['chrom'] == chrom].itertuples(False):
            start = df_chrom.loc[df_chrom['start'] <= row.start].shape[0]
            stop = df_chrom.loc[df_chrom['start'] <= row.stop].shape[0]
            x = range(start, stop + 1)
            ax[i].plot(x, [y] * len(x), color='blue')
            ax[i].text(stop, y, row.name, ha='left', va='center', color='blue')
            y -= 0.2
    plt.savefig(out_file)
    plt.close()
     