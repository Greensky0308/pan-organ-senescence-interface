import os
import numpy as np, pandas as pd

R = 'out'; V = f'{R}/valid'
sen = set(pd.read_csv('raw/sencat_symbols.csv')['symbol'].dropna())
rng = np.random.default_rng(1)

PAIRS = [('kidney', 'kidney'), ('liver', 'liver_cirrhosis'), ('heart', 'heart'), ('lung', 'lung'),
         ('colon', 'colon'), ('skin', 'skin'), ('artery', 'artery'), ('esophagus', 'esophagus'),
         ('brain', 'brain')]

def load(pre):
    e = pd.read_csv(f'{V}/{pre}_expr.tsv', sep='\t', index_col=0)
    g = pd.read_csv(f'{V}/{pre}_groups.tsv', sep='\t').set_index('sample')['group'].reindex(e.columns)
    return e, g

def perm_enrich(bg, sb, up, n=10000):
    obs = len(up & sb)
    if not up: return obs, np.nan, np.nan, 1.0
    perm = np.array([len(set(rng.choice(list(bg), len(up), replace=False)) & sb) for _ in range(n)])
    exp = perm.mean()
    return obs, exp, obs / exp, (np.sum(perm >= obs) + 1) / (n + 1)

rows = []
for org, pre in PAIRS:
    if not (os.path.exists(f'{V}/{pre}_expr.tsv') and os.path.exists(f'{V}/{pre}_up.txt')):
        print(f'{org:9s} MISSING'); continue
    e, _ = load(pre)
    bg = set(e.index.astype(str)); sb = sen & bg
    up = set(x.strip() for x in open(f'{V}/{pre}_up.txt') if x.strip()) & bg
    obs, exp, fold, p = perm_enrich(bg, sb, up)
    rows.append({'organ': org, 'method': 'threshold', 'n_sig': len(up), 'n_bg': len(bg),
                 'overlap': obs, 'exp': round(exp, 2) if exp == exp else np.nan,
                 'fold': round(fold, 2) if fold == fold else np.nan, 'p': p})
    print('%-9s [thr ] sig=%5d | obs=%3d exp=%.1f fold=%.2f p=%.4g' %
          (org, len(up), obs, exp if exp == exp else 0, fold if fold == fold else 0, p))
thr = pd.DataFrame(rows)

print('\n--- top-k rank test (top 5% logFC) ---')
rows2 = []
for org, pre in PAIRS + [('brain28146', 'brain28146')]:
    if not os.path.exists(f'{V}/{pre}_expr.tsv'): continue
    e, g = load(pre)
    bg = set(e.index.astype(str)); sb = sen & bg
    lfc = e.loc[:, g == 'case'].mean(1) - e.loc[:, g == 'control'].mean(1)
    k = max(int(len(bg) * 0.05), 50)
    up = set(lfc.sort_values(ascending=False).head(k).index)
    obs, exp, fold, p = perm_enrich(bg, sb, up, n=5000)
    rows2.append({'organ': org, 'method': 'top5%', 'n_sig': k, 'n_bg': len(bg),
                  'overlap': obs, 'exp': round(exp, 2), 'fold': round(fold, 2), 'p': p})
    print('%-9s [top5%%] k=%4d | obs=%3d exp=%.1f fold=%.2f p=%.4g' % (org, k, obs, exp, fold, p))

out = pd.concat([thr, pd.DataFrame(rows2)], ignore_index=True)
out.to_csv(f'{R}/valid_enrichment.csv', index=False)
print('\nsaved', f'{R}/valid_enrichment.csv')

d = pd.read_csv(f'{R}/interface_enrichment.csv').set_index('organ')
print('\n=== discovery vs validation ===')
for o in d.index:
    t = out[(out.organ == o) & (out.method == 'threshold')]
    if len(t):
        t = t.iloc[0]
        print('%-9s discovery %.2f (%.1g)  ->  validation(threshold) %.2f (%.1g)' % (o, d.loc[o, 'fold'], d.loc[o, 'p'], t['fold'], t['p']))
    else:
        r = out[(out.organ == o) & (out.method == 'top5%')].iloc[0]
        print('%-9s discovery %.2f (%.1g)  ->  validation(top5%%) %.2f (%.1g)' % (o, d.loc[o, 'fold'], d.loc[o, 'p'], r['fold'], r['p']))
