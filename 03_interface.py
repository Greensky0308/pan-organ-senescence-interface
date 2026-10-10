import numpy as np, pandas as pd
from scipy.stats import hypergeom

SEN = 'raw/sencat_symbols.csv'
ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
sen = set(pd.read_csv(SEN)['symbol'].dropna())
print('SenCat unique symbols:', len(sen))

rng = np.random.default_rng(1)
rows, overlap = [], {}
for o in ORGS:
    up = set(x.strip() for x in open(f'out/organ_up_{o}.txt') if x.strip())
    bg = set(pd.read_csv(f'out/cohorts/{o}_expr.tsv', sep='\t', index_col=0, nrows=0).columns)  # placeholder
    bg = set(pd.read_csv(f'out/cohorts/{o}_expr.tsv', sep='\t', usecols=[0], index_col=0).index)
    up_bg = up & bg
    sen_bg = sen & bg
    obs = len(up_bg & sen_bg)
    perm = np.array([len(set(rng.choice(list(bg), len(up_bg), replace=False)) & sen_bg) for _ in range(10000)])
    exp = perm.mean(); p = (np.sum(perm >= obs) + 1) / 10001
    rows.append({'organ': o, 'n_sig': len(up_bg), 'n_bg': len(bg), 'n_sen_in_bg': len(sen_bg),
                 'overlap': obs, 'exp': round(exp, 2), 'fold': round(obs / exp, 2), 'p': p})
    overlap[o] = sorted(up_bg & sen_bg)
    print('%-7s sig=%5d bg=%6d sen_bg=%3d | obs=%3d exp=%.1f fold=%.2f p=%.4g'
          % (o, len(up_bg), len(bg), len(sen_bg), obs, exp, obs / exp, p))

pd.DataFrame(rows).to_csv('out/interface_enrichment.csv', index=False)

allg = sorted(set().union(*overlap.values()))
M = pd.DataFrame(0, index=allg, columns=ORGS)
for o in ORGS:
    M.loc[overlap[o], o] = 1
M['n_organs'] = M.sum(1)
M = M.sort_values(['n_organs'] + ORGS, ascending=False)
M.to_csv('out/interface_matrix.csv')
print('\nTotal interface genes:', len(allg))
print('Counts by number of organs:', M['n_organs'].value_counts().sort_index(ascending=False).to_dict())
print('\nInterface genes in >=4 organs:')
print(M[M['n_organs'] >= 4].to_string())

common_k = {}
for k in range(2, len(ORGS) + 1):
    sets = [set(overlap[o]) for o in ORGS]
    obs_k = sum(1 for g in allg if sum(g in s for s in sets) >= k)
    common_k[k] = obs_k
print('\nNumber of genes shared by >=k organs:', common_k)
