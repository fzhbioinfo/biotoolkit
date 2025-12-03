# -*- coding:utf-8 -*-
import pandas as pd
import glob
import sys
import os
qcs = glob.glob(f'{sys.argv[1]}/F*/*.QC.xlsx')
merge = pd.DataFrame()
for qc in qcs:
    family = os.path.basename(os.path.dirname(qc))
    df = pd.read_excel(qc)
    df['Family'] = family
    merge = merge.append(df, sort=False)
merge.to_excel(sys.argv[2], index=False)
