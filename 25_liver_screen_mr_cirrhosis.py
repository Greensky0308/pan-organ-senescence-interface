import os
import numpy as np, pandas as pd
from scipy.stats import norm

ARC = 'out/causal'
PQTL = {'CHI3L1': 'Folkersen', 'SPP1': 'Folkersen', 'MMP7': 'Folkersen', 'MMP12': 'Folkersen',
        'CCL2': 'Suhre', 'LYZ': 'Suhre', 'TNFAIP6': 'Suhre', 'CXCL1': 'GCST'}
CIR = pd.read_csv(f'{ARC}/cirrhosis_outcome_subset.csv').drop_duplicates('variant_id').set_index('variant_id')

rows = []
for g, src in PQTL.items():
    f = f'{ARC}/coloc_input_scan/{g}__liver.csv'
    if not os.path.exists(f):
        print(g, 'no exposure'); continue
    e = pd.read_csv(f)[['snp', 'ea', 'nea', 'beta.p', 'se.p', 'p']].dropna(subset=['beta.p', 'se.p', 'p']).drop_duplicates('snp')
    gws = e[e['p'] < 5e-8]
    if len(gws) == 0:
        print(g, 'no GWS lead'); continue
    lead = gws.loc[gws['p'].idxmin()]
    if lead['snp'] not in CIR.index:
        print(g, 'lead not in outcome'); continue
    o = CIR.loc[lead['snp']]
    if o['A1'] == lead['ea'] and o['A2'] == lead['nea']:
        by = o['BETA']
    elif o['A1'] == lead['nea'] and o['A2'] == lead['ea']:
        by = -o['BETA']
    else:
        print(g, 'allele mismatch'); continue
    bx = lead['beta.p']; sy = o['SE']
    beta = by / bx; se = abs(sy / bx)
    rows.append({'gene': g, 'source': src, 'organ': 'liver', 'disease': 'Cirrhosis',
                 'lead': lead['snp'], 'lead_p': lead['p'], 'beta': beta,
                 'OR': float(np.exp(beta)), 'p': float(2 * norm.sf(abs(beta / se)))})
R = pd.DataFrame(rows).sort_values('gene')
R.to_csv(f'{ARC}/screen_mr_cirrhosis.csv', index=False)
print(R.to_string(index=False, float_format=lambda x: '%.6g' % x))
