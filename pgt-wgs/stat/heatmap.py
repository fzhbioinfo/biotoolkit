# -*- coding:utf-8 -*-
import sys
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
genome = sys.argv[1]
title = sys.argv[2]
fig = sys.argv[3]
df = pd.read_csv(genome, sep=r'\s+')
df = df[['IID1', 'IID2', 'PI_HAT']].copy()
df = df.pivot('IID1', 'IID2', 'PI_HAT')
#plt.figure(figsize=(20, 18))
ax = sns.heatmap(df, annot=True, vmin=0, vmax=1, cmap="YlGnBu", xticklabels=df.columns.tolist(), yticklabels=df.index.tolist())
plt.title(title)
plt.tight_layout()
plt.savefig(fig)
