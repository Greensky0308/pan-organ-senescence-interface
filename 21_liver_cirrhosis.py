import os, json, time, gzip
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import numpy as np, pandas as pd, requests
from scipy.stats import norm, chi2

ARC = 'out/causal'
RAW = os.path.join(os.environ.get('FINNGEN_DIR', 'data/finngen'), 'finngen_R12_CIRRHOSIS_BROAD.gz')
NCASE, NCTRL = 5545, 494803
N, S = NCASE + NCTRL, NCASE / (NCASE + NCTRL)
LO, HI = 102_300_000, 103_000_000
PICK = {'MMP7': 'Lung', 'MMP12': 'Spleen'}
PQTL = {'MMP7': 'ebi-a-GCST90012056', 'MMP12': 'ebi-a-GCST90012070'}
TOK = open(os.path.expanduser('~/.opengwas_token')).read().strip()
H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}

rows = []
with gzip.open(RAW, 'rt') as fh:
    hdr = fh.readline()
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[0] != '11':
            continue
        pos = int(f[1])
        if pos < LO or pos > HI:
            continue
        rows.append({'snp': f[4].split(',')[0], 'ref': f[2], 'alt': f[3], 'position': pos,
                     'beta': float(f[8]), 'se': float(f[9]), 'eaf': float(f[10])})
OUT = pd.DataFrame(rows).dropna()
OUT['MAF'] = OUT['eaf'].clip(upper=0.5)
OUT['N'] = N; OUT['s'] = S
OUT.to_csv(f'{ARC}/fg_liver_cirrhosis.csv', index=False)
print('cirrhosis cis SNPs:', len(OUT), 'N=%d s=%.4f' % (N, S))

eq = pd.read_csv(f'{ARC}/gtex_eqtl_mmp_independent.csv')
inst = eq[eq.apply(lambda r: PICK.get(r['gene']) == r['tissue'], axis=1)].copy()
inst[['chr', 'pos', 'ref', 'alt']] = inst['variant_id'].str.replace('_b38', '', regex=False).str.split('_', expand=True)
inst['pos'] = inst['pos'].astype(int)
print('\nGTEx instruments:'); print(inst[['gene', 'tissue', 'variant_id', 'slope', 'slope_se', 'pval']].to_string(index=False))
posmap = {(r['position'], r['ref'], r['alt']): (r['beta'], r['se']) for _, r in OUT.iterrows()}
recs = []
for _, s in inst.iterrows():
    key = (s['pos'], s['ref'], s['alt']); swap = (s['pos'], s['alt'], s['ref'])
    if key in posmap:
        by, sy = posmap[key]; flip = 1
    elif swap in posmap:
        by, sy = posmap[swap]; flip = -1
    else:
        continue
    recs.append({'gene': s['gene'], 'bx': s['slope'], 'sx': s['slope_se'], 'by': by * flip, 'sy': sy})
R = pd.DataFrame(recs)
def ivw(gr):
    bx, sx, by, sy = gr['bx'].values, gr['sx'].values, gr['by'].values, gr['sy'].values
    ratios = by / bx; w = bx**2 / sy**2
    beta = np.sum(w * ratios) / np.sum(w); se = np.sqrt(1 / np.sum(w))
    Q = np.sum(w * (ratios - beta)**2)
    return pd.Series({'n_snp': len(gr), 'beta': beta, 'se': se, 'p': 2 * norm.sf(abs(beta / se)),
                      'Q': Q, 'Q_p': chi2.sf(Q, len(gr) - 1) if len(gr) > 1 else np.nan})
E = R.groupby('gene').apply(ivw, include_groups=False).reset_index()
E['organ'] = 'liver'; E['disease'] = 'Cirrhosis'; E['OR'] = np.exp(E['beta'])
E = E[['organ', 'disease', 'gene', 'n_snp', 'beta', 'se', 'p', 'Q', 'Q_p', 'OR']]
E.to_csv(f'{ARC}/mr_eqtl_cirrhosis.csv', index=False)
print('\n=== MR eQTL -> cirrhosis ==='); print(E.to_string(index=False, float_format=lambda x: '%.3g' % x))

sess = requests.Session(); sess.proxies = ({'https': os.environ['HTTPS_PROXY']} if os.environ.get('HTTPS_PROXY') else None)
def get_pqtl(dsid):
    for _ in range(6):
        try:
            j = sess.post('https://api.opengwas.io/api/associations', headers=H,
                          data=json.dumps({'id': [dsid], 'variant': ['11:102300000-103000000']}), timeout=120).json()
            return pd.DataFrame(j)
        except Exception:
            time.sleep(2)
    raise RuntimeError('pQTL fetch failed ' + dsid)

def clump(df, window=1_000_000):
    df = df.sort_values('p'); keep = []
    for _, r in df.iterrows():
        if all(abs(r['position'] - k) > window for k in keep):
            keep.append(r['position'])
    return df[df['position'].isin(keep)]

orec = []
for g, ds in PQTL.items():
    pq = get_pqtl(ds)
    pq = pq.dropna(subset=['rsid', 'beta', 'se', 'eaf', 'p'])
    pq = pq.rename(columns={'rsid': 'snp'})
    pq['MAF'] = pq['eaf'].clip(upper=0.5)
    pq[['snp', 'beta', 'se', 'eaf', 'MAF', 'n', 'ea', 'nea', 'p', 'position']].to_csv(f'{ARC}/pqtl_{g}_cirrhosis.csv', index=False)
    sig = pq[pq['p'] < 5e-8].copy(); sig = clump(sig)
    print('\npQTL %-5s: cis SNP=%d, GWS=%d, after clump=%d (lead %s p=%.2g)'
          % (g, len(pq), (pq['p'] < 5e-8).sum(), len(sig), sig.iloc[0]['snp'], sig.iloc[0]['p']))
    omap = OUT.set_index('snp')
    for _, s in sig.iterrows():
        if s['snp'] not in omap.index:
            continue
        o = omap.loc[s['snp']]
        if o['alt'] == s['ea'] and o['ref'] == s['nea']:
            b = o['beta']
        elif o['alt'] == s['nea'] and o['ref'] == s['ea']:
            b = -o['beta']
        else:
            continue
        orec.append({'organ': 'liver', 'disease': 'Cirrhosis', 'gene': g, 'rsid': s['snp'],
                     'bx': s['beta'], 'sx': s['se'], 'by': b, 'sy': o['se']})
P = pd.DataFrame(orec)
Pp = P.groupby(['gene', 'organ', 'disease']).apply(ivw, include_groups=False).reset_index()
Pp['OR'] = np.exp(Pp['beta'])
Pp = Pp[['gene', 'organ', 'disease', 'n_snp', 'beta', 'se', 'p', 'Q', 'Q_p', 'OR']]
Pp.to_csv(f'{ARC}/mr_pqtl_cirrhosis.csv', index=False)
print('\n=== MR pQTL -> cirrhosis ==='); print(Pp.to_string(index=False, float_format=lambda x: '%.3g' % x))
