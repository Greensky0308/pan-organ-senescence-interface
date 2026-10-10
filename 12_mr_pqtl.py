import os, time, json
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import numpy as np, pandas as pd, requests, pysam
from scipy.stats import norm, chi2

TOK = os.environ.get('OPENGWAS_JWT', '')
H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}
PQTL = {'MMP7': 'ebi-a-GCST90012056', 'MMP12': 'ebi-a-GCST90012070'}   # Folkersen plasma pQTL
CIS = '11:102300000-103000000'                                        # hg19 cis window
DISEASE = {'kidney': ('DKD','DM_NEPHROPATHY'), 'liver': ('MASLD','NAFLD'),
           'heart': ('DCM','I9_NONISCHCARDMYOP'), 'lung': ('IPF','IPF'),
           'brain': ('AD','G6_ALZHEIMER'), 'colon': ('IBD','K11_IBD_STRICT'),
           'skin': ('Psoriasis','L12_PSORIASIS'), 'artery': ('Athero','I9_ATHSCLE'),
           'esophagus': ('EAC','C3_OESOPHAGUS_EXALLC')}
FF = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz'

def get_pqtl(dsid):
    s = requests.Session(); s.proxies = ({'https': os.environ['HTTPS_PROXY']} if os.environ.get('HTTPS_PROXY') else None)
    for _ in range(5):
        try:
            r = s.post('https://api.opengwas.io/api/associations', headers=H,
                       data=json.dumps({'id': [dsid], 'variant': [CIS]}), timeout=120)
            return pd.DataFrame(r.json())
        except Exception:
            time.sleep(2)
    raise RuntimeError('pQTL fetch failed')

def fetch_lines(url, contig, lo, hi, tries=8):
    last = None
    for _ in range(tries):
        try:
            return list(pysam.TabixFile(url).fetch(contig, lo, hi))
        except Exception as e:
            last = e; time.sleep(1.5)
    raise last

def clump(df, window=1_000_000):
    df = df.sort_values('p')
    keep = []
    for _, r in df.iterrows():
        if all(abs(r['position'] - k) > window for k in keep):
            keep.append(r['position'])
    return df[df['position'].isin(keep)]

inst = {}
for g, ds in PQTL.items():
    pq = get_pqtl(ds)
    sig = pq[pq['p'] < 5e-8].copy()
    sig = clump(sig)
    inst[g] = sig
    print('%-6s %s: cis SNP=%d, GWS=%d, after clump=%d (lead rsid %s p=%.2g)'
          % (g, ds, len(pq), (pq['p'] < 5e-8).sum(), len(sig), sig.iloc[0]['rsid'], sig.iloc[0]['p']))

def fg_index(dsid_pheno):
    fg = {}
    for line in fetch_lines(FF % dsid_pheno, '11', 102_300_000, 103_000_000):
        f = line.split('\t')  # chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta ...
        for rs in f[4].split(','):
            fg[rs] = (f[2], f[3], float(f[8]), float(f[9]))  # ref, alt, beta(per alt), sebeta
    return fg

recs = []
for organ, (lab, pheno) in DISEASE.items():
    fg = fg_index(pheno)
    for g, sig in inst.items():
        for _, s in sig.iterrows():
            if s['rsid'] not in fg: continue
            ref, alt, by, sy = fg[s['rsid']]
            if alt == s['ea'] and ref == s['nea']:   b = by
            elif alt == s['nea'] and ref == s['ea']: b = -by
            else: continue
            recs.append({'organ': organ, 'disease': lab, 'gene': g, 'rsid': s['rsid'],
                         'bx': s['beta'], 'sx': s['se'], 'by': b, 'sy': sy})
R = pd.DataFrame(recs)

def mr(gr):
    bx, sx, by, sy = gr['bx'].values, gr['sx'].values, gr['by'].values, gr['sy'].values
    ratios = by / bx; w = bx**2 / sy**2
    beta = np.sum(w * ratios) / np.sum(w); se = np.sqrt(1 / np.sum(w))
    p = 2 * norm.sf(abs(beta / se))
    Q = np.sum(w * (ratios - beta)**2)
    return pd.Series({'n_snp': len(gr), 'beta': beta, 'se': se, 'p': p,
                      'Q': Q, 'Q_p': chi2.sf(Q, len(gr)-1) if len(gr) > 1 else np.nan})

res = R.groupby(['gene', 'organ', 'disease']).apply(mr, include_groups=False).reset_index()
res['OR'] = np.exp(res['beta'])
res = res.sort_values(['gene', 'p'])
res.to_csv('out/mr_pqtl_results.csv', index=False)
print('\n=== pQTL MR (Folkersen → FinnGen) ===')
print(res.to_string(index=False, float_format=lambda x: '%.3g' % x))
