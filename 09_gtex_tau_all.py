import numpy as np, pandas as pd

hg = pd.read_csv('raw/hgnc.txt', sep='\t', dtype=str, low_memory=False)
ENS2SYM = dict(zip(hg['ensembl_gene_id'].str.split('.').str[0], hg['symbol']))
g = pd.read_csv('raw/gtex_median.gct.gz', sep='\t', skiprows=2)
g['sym'] = [ENS2SYM.get(str(n).split('.')[0]) for n in g['Name']]
g = g[g['sym'].notna()].groupby('sym').mean(numeric_only=True)
tpm = g.values
mx = np.maximum(tpm.max(1, keepdims=True), 1e-9)
tau = (1 - tpm / mx).sum(1) / (tpm.shape[1] - 1)

M = pd.read_csv('out/interface_matrix.csv', index_col=0)
iface = set(M.index)
out = pd.DataFrame({'tau': tau, 'is_iface': [ix in iface for ix in g.index]}, index=g.index)
out.index.name = 'gene'
out = out.dropna(subset=['tau'])
out.to_csv('out/gtex_tau_all.csv')
print(out.groupby('is_iface')['tau'].agg(['size', 'median']))
