# -*- coding:utf-8 -*-
from bioutils.assemblies import make_ac_name_map
from io import StringIO
import pandas as pd
import numpy as np
import subprocess
import sys


def ac_to_name(chromosome_dic, chromosome):
    name = chromosome_dic[chromosome]
    if not name.startswith('chr'):
        name = 'chr' + name
    return name


gff = sys.argv[1]
genome = sys.argv[2]
chrome_dic = make_ac_name_map(genome)
chrome_dic['NC_012920.1'] = 'chrM_NC_012920.1'
# read and extract info
annotation = subprocess.getoutput('zcat ' + gff + ' | grep -v "#"')
annotation_df = pd.read_csv(StringIO(annotation), sep='\t', header=None)
annotation_df['gene'] = annotation_df[8].str.extract('gene=(.*?);')
annotation_df['ID'] = annotation_df[8].str.extract('ID=(.*?);')
annotation_df['tag'] = annotation_df[8].str.extract('tag=(.*?);')
annotation_df['transcript_id'] = annotation_df[8].str.extract('transcript_id=(.*?)$')
annotation_df.fillna('.', inplace=True)
# filter RefSeq Select exon info
annotation_df = annotation_df[annotation_df[2].isin(['exon']) & annotation_df[1].str.contains('RefSeq') & annotation_df[0].str.contains('NC_')].copy()
annotation_df_select = annotation_df[annotation_df['tag'] == 'RefSeq Select'].copy()
# process Not RefSeq Select, choose longest
select_genes = annotation_df_select['gene'].unique().tolist()
annotation_df_notselect = annotation_df[(~annotation_df['gene'].isin(select_genes)) & (annotation_df['transcript_id'] != '.')].copy()
notselect = annotation_df_notselect.copy()
notselect['size'] = notselect[4] - notselect[3] + 1
notselect_stat = notselect.groupby(['gene', 'transcript_id']).agg({'size': np.sum}).reset_index()
notselect_stat.sort_values(by='size', ascending=False, inplace=True)
notselect_stat.drop_duplicates(subset=['gene'], keep='first', inplace=True)
annotation_df_notselect_choose = annotation_df_notselect[annotation_df_notselect['transcript_id'].isin(notselect_stat['transcript_id'].unique().tolist())].copy()
# final annotation info
annotation_df_final = annotation_df_select.append(annotation_df_notselect_choose)
annotation_df_final[0] = annotation_df_final.apply(lambda x: ac_to_name(chrome_dic, x[0]), axis=1)
annotation_df_final.to_csv(sys.argv[3], sep='\t', columns=[0, 3, 4, 6, 'gene', 'ID', 'tag'], index=False, header=None)