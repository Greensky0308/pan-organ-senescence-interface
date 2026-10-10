import numpy as np, pandas as pd
from collections import Counter
import os

ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
CFG = {'kidney': (['DN', 'FSGS', 'MCD', 'RPGN'], 2), 'colon': (['ulcerative_colitis', 'CD'], 1),
       'skin': (['case'], 1), 'lung': (['ipf'], 1),
       'liver': (['cirrhosis'], 1),
       'heart': (['Dilated_cardiomyopathy__DCM_'], 1), 'artery': (['case'], 1), 'esophagus': (['case'], 1), 'brain': (['case'], 1)}

def load_sig(fn):
    return set(x.strip() for x in open(fn) if x.strip() and len(x.strip()) > 1)

SIGS = {'SenCat': set(pd.read_csv('raw/sencat_symbols.csv')['symbol'].dropna())}
SIGS['Reactome'] = load_sig('raw/sig_reactome.txt')
SIGS['GO0090398'] = load_sig('raw/sig_go_senescence.txt')
SIGS['SenMayo'] = load_sig('raw/sig_senmayo.txt')
SIGS['Fridman'] = load_sig('raw/sig_fridman.txt')
SIGS['ReactomeSASP'] = load_sig('raw/sig_reactome_sasp.txt')
for k, v in SIGS.items():
    print(k, len(v))

ups, bgs = {}, {}
for o in ORGS:
    bg = set(pd.read_csv(f'out/cohorts/{o}_expr.tsv', sep='\t', usecols=[0], index_col=0).index)
    cs, k = CFG[o]
    cnt = Counter()
    for c in cs:
        f = f'out/cohorts/{o}_DE_{c}.csv'
        if not os.path.exists(f): continue
        d = pd.read_csv(f)
        for g in d.loc[(d.FDR < 0.05) & (d.logFC >= 1), 'gene']:
            cnt[g] += 1
    ups[o] = {g for g, n in cnt.items() if n >= k} & bg
    bgs[o] = bg

rng = np.random.default_rng(2)
rows = []
for sname, S in SIGS.items():
    for o in ORGS:
        sb = S & bgs[o]; u = ups[o] & bgs[o]
        obs = len(u & sb)
        perm = np.array([len(set(rng.choice(list(bgs[o]), len(u), replace=False)) & sb) for _ in range(4000)])
        exp = perm.mean(); p = (np.sum(perm >= obs) + 1) / 4001
        rows.append({'sig': sname, 'organ': o, 'n_sig': len(u), 'n_sen_bg': len(sb),
                     'obs': obs, 'exp': round(exp, 2), 'fold': round(obs / exp, 2) if exp > 0 else np.nan, 'p': p})
df = pd.DataFrame(rows)
df.to_csv('out/multisig_enrichment.csv', index=False)
print('\n=== Enrichment fold per signature ===')
piv = df.pivot(index='organ', columns='sig', values='fold').loc[ORGS]
print(piv.round(2).to_string())
piv_p = df.pivot(index='organ', columns='sig', values='p').loc[ORGS]
print('\nSignificant (p<0.05):', {s: int((piv_p[s] < 0.05).sum()) for s in piv_p.columns})

print('\n=== Core size (>=3 organs) per signature ===')
for sname, S in SIGS.items():
    ov = {o: ups[o] & bgs[o] & S for o in ORGS}
    allg = set().union(*ov.values())
    n3 = sum(1 for g in allg if sum(g in ov[o] for o in ORGS) >= 3)
    n5 = sum(1 for g in allg if sum(g in ov[o] for o in ORGS) >= 5)
    mmp7 = sum(1 for o in ORGS if 'MMP7' in ov[o])
    print(f'{sname:10s} interface {len(allg):4d} | >=3 organs {n3:3d} | >=5 organs {n5} | MMP7 in {mmp7} organs')
