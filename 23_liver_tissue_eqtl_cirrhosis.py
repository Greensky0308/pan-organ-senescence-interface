import gzip, os
import numpy as np, pandas as pd
from scipy.stats import norm

ARC = 'out/causal'
RAW = os.path.join(os.environ.get('FINNGEN_DIR', 'data/finngen'), 'finngen_R12_CIRRHOSIS_BROAD.gz')

inst = pd.read_csv(f'{ARC}/gtex_core_tissue_eqtl.csv')
inst = inst[inst['tissue'] == 'Liver'].copy()
inst[['chr', 'pos', 'ref', 'alt']] = inst['variant_id'].str.replace('_b38', '', regex=False).str.split('_', expand=True)
inst['pos'] = inst['pos'].astype(int)
inst['chr'] = inst['chr'].str.replace('chr', '', regex=False)

need = {}
for _, s in inst.iterrows():
    need.setdefault(s['chr'], []).append(s['pos'])
ranges = {c: (min(v) - 150_000, max(v) + 150_000) for c, v in need.items()}
print('windows:', ranges)

rows = []
with gzip.open(RAW, 'rt') as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip('\n').split('\t')
        c = f[0]
        if c not in ranges:
            continue
        pos = int(f[1])
        lo, hi = ranges[c]
        if pos < lo or pos > hi:
            continue
        rows.append({'snp': f[4].split(',')[0], 'ref': f[2], 'alt': f[3], 'position': pos,
                     'beta': float(f[8]), 'se': float(f[9])})
OUT = pd.DataFrame(rows)
pmap = {(r['position'], r['ref'], r['alt']): (r['beta'], r['se']) for _, r in OUT.iterrows()}

recs = []
for _, s in inst.iterrows():
    key = (s['pos'], s['ref'], s['alt']); swap = (s['pos'], s['alt'], s['ref'])
    if key in pmap:
        by, sy = pmap[key]; flip = 1
    elif swap in pmap:
        by, sy = pmap[swap]; flip = -1
    else:
        print('miss', s['gene'], s['variant_id']); continue
    bx = s['slope']; by *= flip
    beta = by / bx; se = sy / abs(bx)
    recs.append({'gene': s['gene'], 'organ': 'liver', 'disease': 'Cirrhosis', 'tissue': 'Liver',
                 'snp': s['variant_id'], 'eqtl_p': s['pval'], 'beta': beta,
                 'OR': float(np.exp(beta)), 'p': 2 * norm.sf(abs(beta / se))})
R = pd.DataFrame(recs)
R.to_csv(f'{ARC}/eqtl_tissue_mr_cirrhosis.csv', index=False)
print(R.to_string(index=False, float_format=lambda x: '%.6g' % x))
