from bioutils.assemblies import make_ac_name_map
import pandas as pd
import sys


def return_ac_name(ac_name_dic, ac):
    return f'chr{ac_name_dic[ac]}'


gap = pd.read_csv(sys.argv[1], sep='\t')
version = sys.argv[2]
dic = make_ac_name_map(version)
gap = gap[gap['# accession.version'].str.startswith('NC_')].copy()
gap['# accession.version'] = gap.apply(lambda x: return_ac_name(dic, x['# accession.version']), axis=1)
gap.to_csv(sys.argv[3], sep='\t', index=False)
