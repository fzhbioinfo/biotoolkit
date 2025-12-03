from collections import defaultdict
from pyfaidx import Fasta
import pandas as pd
import numpy as np
import sys
import re


def get_n_size(fa, chromosome, start, stop):
    if start == 0:
        start = 1
    return len(re.findall(r'[nN]', str(fa.get_seq(chromosome, start, stop))))


coverage = pd.read_csv(sys.argv[1], sep='\t')
ref = Fasta(sys.argv[2])
coverage['N_size'] = coverage.apply(lambda x: get_n_size(ref, x['chromosome'], x['start'], x['stop']), axis=1)
print(coverage[['chromosome', 'start', 'stop', 'size', 'N_size']])
coverage.drop(columns=['chromosome', 'start', 'stop'], inplace=True)
dic = coverage.sum().to_dict()
stat = defaultdict(dict)
stat['mean_depth(rm_dup)'] = np.round(dic['base_on_target(rm_dup)'] / (dic['size'] - dic['N_size']), 4)
stat['coverage>=1(rm_dup)'] = np.round(dic['coverage>=1(rm_dup)'] / (dic['size'] - dic['N_size']), 4)
stat['coverage>=4(rm_dup)'] = np.round(dic['coverage>=4(rm_dup)'] / (dic['size'] - dic['N_size']), 4)
stat['coverage>=10(rm_dup)'] = np.round(dic['coverage>=10(rm_dup)'] / (dic['size'] - dic['N_size']), 4)
stat['coverage>=20(rm_dup)'] = np.round(dic['coverage>=20(rm_dup)'] / (dic['size'] - dic['N_size']), 4)
stat_df = pd.DataFrame.from_dict(stat, orient='index').T
stat_df.to_csv(sys.argv[3], sep='\t', index=False)
