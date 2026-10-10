import gzip, random
import numpy as np, pandas as pd
from scipy.stats import spearmanr

AGES = {'20-29': 1, '30-39': 2, '40-49': 3, '50-59': 4, '60-69': 5, '70-79': 6}
meta = pd.read_csv('out/gtex_samples.csv')
meta = meta[meta['age'].isin(AGES)]
meta['a'] = meta['age'].map(AGES)
samp2tis = dict(zip(meta['sample'], meta['tissue']))
samp2age = dict(zip(meta['sample'], meta['a']))

M = pd.read_csv('out/interface_matrix.csv', index_col=0)
IFACE = set(M.index)

# HGNC ens->sym
hg = pd.read_csv('raw/hgnc.txt', sep='\t', dtype=str, low_memory=False)
ENS2SYM = dict(zip(hg['ensembl_gene_id'].str.split('.').str[0], hg['symbol']))

random.seed(0)
keep, rows = {}, {}
want_bg = set()
with gzip.open('raw/gtex_tpm.gct.gz', 'rt') as fh:
    fh.readline(); fh.readline()
    hdr = fh.readline().rstrip('\n').split('\t')
    cols = hdr[2:]
    rng = random.Random(1)
    for line in fh:
        p = line.rstrip('\n').split('\t')
        sym = ENS2SYM.get(p[0].split('.')[0])
        if sym is None: continue
        if sym in IFACE:
            rows[sym] = np.array(p[2:], dtype='float32')
        elif sym not in keep and rng.random() < 0.06:
            keep[sym] = 1; rows[sym] = np.array(p[2:], dtype='float32')
print('kept genes:', len(rows), '| interface:', sum(1 for s in rows if s in IFACE))

expr = pd.DataFrame(rows, index=cols).T          # gene x sample
tis_samples = {}
for s in cols:
    if s in samp2tis: tis_samples.setdefault(samp2tis[s], []).append(s)
res = []
for tis, ss in tis_samples.items():
    if len(ss) < 40: continue
    age = np.array([samp2age[s] for s in ss])
    if len(set(age)) < 3: continue
    E = expr[ss]
    for g in E.index:
        x = E.loc[g].values
        if np.nanstd(x) == 0: continue
        rho, p = spearmanr(x, age, nan_policy='omit')
        if not np.isnan(rho):
            res.append({'gene': g, 'tissue': tis, 'rho': rho, 'p': p, 'n': len(ss)})
R = pd.DataFrame(res)
R.to_csv('out/gtex_age_cor.csv', index=False)
summ = R.groupby('gene').agg(mean_rho=('rho', 'mean'), n_pos=('rho', lambda v: (v > 0).sum()),
                             n_tis=('rho', 'size'), n_sig_pos=('p', lambda v: 0)).reset_index()
sig = R[(R['p'] < 0.05)]
pos_sig = sig[sig['rho'] > 0].groupby('gene').size()
summ['n_sig_pos'] = summ['gene'].map(pos_sig).fillna(0).astype(int)
summ['is_iface'] = summ['gene'].isin(IFACE)
summ.to_csv('out/gtex_age_summary.csv', index=False)
print('\nInterface vs background age association:')
print(summ.groupby('is_iface').agg(mean_rho=('mean_rho', 'mean'), n_sig_pos=('n_sig_pos', 'mean'),
                                   n=('gene', 'size')).round(3).to_string())
core = list(M[M['n_organs'] >= 4].index)
print('\nCore (>=4 organs) age association:')
print(summ[summ['gene'].isin(core)][['gene', 'mean_rho', 'n_sig_pos', 'n_tis']].to_string(index=False))
