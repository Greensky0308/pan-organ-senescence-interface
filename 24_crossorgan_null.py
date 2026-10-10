import json
import numpy as np, pandas as pd

ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
M = pd.read_csv('out/interface_matrix.csv', index_col=0)
iface = {o: set(M.index[M[o] == 1]) for o in ORGS}
bg = {}
for o in ORGS:
    cols = pd.read_csv(f'out/cohorts/{o}_expr.tsv', sep='\t', usecols=[0], index_col=0).index
    bg[o] = set(cols)

ks = [2, 3, 4, 5, 6]
obs = {k: int((M['n_organs'] >= k).sum()) for k in ks}
rng = np.random.default_rng(1)
N = 1000
null = {k: [] for k in ks}
for _ in range(N):
    samp = []
    for o in ORGS:
        b = list(bg[o])
        n = len(iface[o])
        if n > len(b): n = len(b)
        samp.append(set(rng.choice(b, size=n, replace=False)))
    from collections import Counter
    c = Counter()
    for s in samp: c.update(s)
    for k in ks:
        null[k].append(sum(1 for g, v in c.items() if v >= k))
out = {'obs': {str(k): obs[k] for k in ks},
       'nullmean': {str(k): float(np.mean(null[k])) for k in ks},
       'p': {str(k): float((1 + np.sum(np.array(null[k]) >= obs[k])) / (N + 1)) for k in ks}}
json.dump(out, open('out/crossorgan_null.json', 'w'))
print('obs   ', out['obs'])
print('null  ', {k: round(v, 2) for k, v in out['nullmean'].items()})
print('p     ', {k: round(v, 4) for k, v in out['p'].items()})
