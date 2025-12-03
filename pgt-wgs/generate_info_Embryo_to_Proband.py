import pandas as pd
import sys
info = pd.read_csv(sys.argv[1], sep='\t')
proband = info.loc[info['sample'].str.contains('Embryo')].copy()
proband['sample'] = 'Proband'
proband['sampleID'] = proband['sampleID'] + '-p'
proband['relation'] = '先证者'
info.append(proband).to_csv(sys.argv[2], sep='\t', index=False)

