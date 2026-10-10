import os
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import re, time
import numpy as np, pandas as pd
import pysam
from scipy.stats import norm, chi2

def fetch_lines(url, contig, lo, hi, tries=8):
    last = None
    for _ in range(tries):
        try:
            return list(pysam.TabixFile(url).fetch(contig, lo, hi))
        except Exception as e:
            last = e; time.sleep(1.5)
    raise last

BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz'
DISEASE = {
    'kidney':    ('DKD',       'DM_NEPHROPATHY'),
    'liver':     ('MASLD',     'NAFLD'),
    'heart':     ('DCM',       'I9_NONISCHCARDMYOP'),
    'lung':      ('IPF',       'IPF'),
    'brain':     ('AD',        'G6_ALZHEIMER'),
    'colon':     ('IBD',       'K11_IBD_STRICT'),
    'skin':      ('Psoriasis', 'L12_PSORIASIS'),
    'artery':    ('Athero',    'I9_ATHSCLE'),
    'esophagus': ('EAC',       'C3_OESOPHAGUS_EXALLC'),
}
PICK = {'MMP7': 'Lung', 'MMP12': 'Spleen'}

eq = pd.read_csv('out/gtex_eqtl_mmp_independent.csv')
inst = eq[eq.apply(lambda r: PICK.get(r['gene']) == r['tissue'], axis=1)].copy()
inst[['chr', 'pos', 'ref', 'alt']] = inst['variant_id'].str.replace('_b38', '', regex=False).str.split('_', expand=True)
inst['pos'] = inst['pos'].astype(int)
print('instrument SNPs:'); print(inst[['gene', 'tissue', 'variant_id', 'slope', 'slope_se', 'pval']].to_string(index=False))

LO, HI = 102_500_000, 103_000_000
recs = []
for organ, (lab, pheno) in DISEASE.items():
    tab = BASE % pheno
    fg = {}
    for line in fetch_lines(tab, '11', LO, HI):
        f = line.split('\t')
        fg[(int(f[1]), f[2], f[3])] = (float(f[8]), float(f[9]))  # beta, sebeta keyed by pos/ref/alt
    for _, s in inst.iterrows():
        key = (s['pos'], s['ref'], s['alt'])
        swap = (s['pos'], s['alt'], s['ref'])
        if key in fg:
            by, sy = fg[key]; flip = 1
        elif swap in fg:
            by, sy = fg[swap]; flip = -1
        else:
            continue
        recs.append({'organ': organ, 'disease': lab, 'gene': s['gene'],
                     'bx': s['slope'], 'sx': s['slope_se'], 'by': by * flip, 'sy': sy,
                     'p_exp': s['pval'], 'p_out': None})
R = pd.DataFrame(recs)

def mr(gr):
    bx, sx, by, sy = gr['bx'].values, gr['sx'].values, gr['by'].values, gr['sy'].values
    ratios = by / bx
    w = bx**2 / sy**2
    beta = np.sum(w * ratios) / np.sum(w)
    se = np.sqrt(1 / np.sum(w))
    z = beta / se; p = 2 * norm.sf(abs(z))
    Q = np.sum(w * (ratios - beta)**2)    # Cochran's Q
    Qp = chi2.sf(Q, len(gr) - 1) if len(gr) > 1 else np.nan
    return pd.Series({'n_snp': len(gr), 'beta': beta, 'se': se, 'p': p, 'Q': Q, 'Q_p': Qp})

res = R.groupby(['organ', 'disease', 'gene']).apply(mr, include_groups=False).reset_index()
res['OR'] = np.exp(res['beta'])
res = res.sort_values(['gene', 'p'])
res.to_csv('out/mr_eqtl_results.csv', index=False)
pd.set_option('display.width', 200)
print('\n=== MR (eQTL → FinnGen) ===')
print(res.to_string(index=False, float_format=lambda x: '%.3g' % x))
