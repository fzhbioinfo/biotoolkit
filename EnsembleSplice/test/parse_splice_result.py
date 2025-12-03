# -*- coding:utf-8 -*-
from io import StringIO
import pandas as pd
import numpy as np
import subprocess
import sys


def parse_maxentscan(result):
    try:
        return -float(result.split('|')[-1])
    except:
        return None


def parse_scsnv(result):
    ada_score = result.split('|')[-2]
    rf_score = result.split('|')[-1]
    if '.' not in [ada_score, rf_score]:
        return np.mean([float(ada_score), float(rf_score)])
    elif ada_score == '.' and rf_score != '.':
        return float(rf_score)
    elif ada_score != '.' and rf_score == '.':
        return float(ada_score)
    else:
        return None


def parse_spliceai(result):
    score = max(result.split('|')[2:6])
    if score == '.':
        return None
    return score


vcf = pd.read_csv(StringIO(subprocess.getoutput("grep -v '##' " + sys.argv[1])), sep='\t')
vcf.rename(columns={'ID': 'Effect'}, inplace=True)
prediction = vcf['INFO'].str.split(';', expand=True)
prediction['MaxEntScan'] = prediction.apply(lambda x: parse_maxentscan(x[0]), axis=True)
prediction['dbscSNV'] = prediction.apply(lambda x: parse_scsnv(x[1]), axis=True)
prediction['SpliceAI'] = prediction.apply(lambda x: parse_spliceai(x[2]), axis=True)
df = vcf.join(prediction)
df.to_csv(sys.argv[2], sep='\t', index=False, columns=['#CHROM', 'POS', 'REF', 'ALT', 'Effect', 'MaxEntScan', 'dbscSNV', 'SpliceAI'])
