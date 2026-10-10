import os
import numpy as np, pandas as pd

ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
SUBTYPE = {
    'kidney':    ['DN', 'FSGS', 'MCD', 'RPGN'],
    'liver':     ['cirrhosis'],
    'colon':     ['CD', 'ulcerative_colitis'],
    'heart':     ['Dilated_cardiomyopathy__DCM_'],
    'lung':      ['ipf'],
    'brain':     ['case'],
    'skin':      ['case'],
    'artery':    ['case'],
    'esophagus': ['case'],
}
UP = lambda df: set(df.loc[(df['FDR'] < 0.05) & (df['logFC'] >= 1), 'gene'])

M = pd.read_csv('out/interface_matrix.csv', index_col=0)

rows = []
ratios_iface, ratios_bg = [], []
per_gene = []
for o in ORGS:
    subs = []
    universe = None
    for s in SUBTYPE[o]:
        f = f'out/cohorts/{o}_DE_{s}.csv'
        df = pd.read_csv(f)
        subs.append(UP(df))
        universe = set(df['gene']) if universe is None else (universe & set(df['gene']))
    iface = [g for g in M.index[M[o] == 1] if g in universe]
    bg = list(universe - set(M.index))
    n_sub = len(subs)

    def ratio(vec):
        return np.array([sum(g in u for u in subs) / n_sub for g in vec])

    ri = ratio(iface); rb = ratio(bg)
    ratios_iface.append(ri); ratios_bg.append(rb)
    for g, r in zip(iface, ri):
        per_gene.append({'organ': o, 'gene': g, 'n_subtype': n_sub, 'n_up': int(round(r * n_sub)), 'ratio': r})
    rows.append({'organ': o, 'n_subtype': n_sub, 'n_iface': len(iface),
                 'iface_ratio_med': float(np.median(ri)), 'iface_ratio_mean': float(ri.mean()),
                 'bg_ratio_med': float(np.median(rb)), 'bg_ratio_mean': float(rb.mean())})
    print('%-9s subs=%d iface(n=%3d) med=%.2f mean=%.2f | bg(n=%5d) med=%.2f mean=%.2f'
          % (o, n_sub, len(iface), np.median(ri), ri.mean(), len(bg), np.median(rb), rb.mean()))

pd.DataFrame(rows).to_csv('out/subtype_consistency_summary.csv', index=False)
PG = pd.DataFrame(per_gene)
PG.to_csv('out/subtype_consistency_per_gene.csv', index=False)

multi = PG[PG['n_subtype'] > 1]
print('\n--- Interface-gene consistency across subtypes ---')
for o in ['kidney', 'colon']:
    s = multi[multi['organ'] == o]
    print('%-9s n=%3d  up ratio med=%.2f  100%% consistent=%d (%.0f%%)'
          % (o, len(s), s['ratio'].median(), (s['ratio'] == 1).sum(), 100 * (s['ratio'] == 1).mean()))

core = list(M[M['n_organs'] >= 4].index)
print('\n--- Core-gene subtype coverage ---')
for g in core:
    parts = []
    for o in ORGS:
        if M.loc[g, o] == 1:
            r = multi[(multi['gene'] == g) & (multi['organ'] == o)]
            if len(r):
                parts.append('%s %d/%d' % (o, r['n_up'].iloc[0], r['n_subtype'].iloc[0]))
            else:
                parts.append('%s 1/1' % o)
    print('%-8s %s' % (g, ' | '.join(parts)))
