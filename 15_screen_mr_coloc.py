import os, json, time
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import numpy as np, pandas as pd, requests, pysam
from scipy.stats import norm

TOK = os.environ.get('OPENGWAS_JWT', '')
H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}
S = requests.Session(); S.proxies = ({'https': os.environ['HTTPS_PROXY']} if os.environ.get('HTTPS_PROXY') else None)
OUT = 'out/coloc_input_scan'; CACHE = 'cache/fgcache'; os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)

COORD = json.load(open('cache/gene_coords.json'))
PQTL = {'CHI3L1': ('ebi-a-GCST90012039', 21758, 'Folkersen'),
        'SPP1':   ('ebi-a-GCST90010244', 21758, 'Folkersen'),
        'MMP7':   ('ebi-a-GCST90012056', 21758, 'Folkersen'),
        'MMP12':  ('ebi-a-GCST90012070', 21758, 'Folkersen'),
        'CCL2':   ('prot-b-4', 3394, 'Suhre'),
        'LYZ':    ('prot-a-1821', 3301, 'Suhre'),
        'TNFAIP6':('prot-a-3031', 3301, 'Suhre'),
        'CXCL1':  ('ebi-a-GCST004457', 3505, 'GCST')}
DISEASE = {'kidney':('DKD','DM_NEPHROPATHY'),'liver':('MASLD','NAFLD'),'heart':('DCM','I9_NONISCHCARDMYOP'),
           'lung':('IPF','IPF'),'brain':('AD','G6_ALZHEIMER'),'colon':('IBD','K11_IBD_STRICT'),
           'skin':('Psoriasis','L12_PSORIASIS'),'artery':('Athero','I9_ATHSCLE'),'esophagus':('EAC','C3_OESOPHAGUS_EXALLC')}
FF = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz'
WIN = 500_000

def og_cis(dsid, chrom, lo, hi, tries=8):
    var = '%s:%d-%d' % (chrom, int(lo), int(hi))
    for _ in range(tries):
        try:
            j = S.post('https://api.opengwas.io/api/associations', headers=H,
                       data=json.dumps({'id': [dsid], 'variant': [var]}), timeout=120).json()
            if isinstance(j, list) and j and 'beta' in j[0]:
                return pd.DataFrame(j)
        except Exception:
            pass
        time.sleep(3)
    return pd.DataFrame()

def fg_cis(pheno, gene, chrom, lo, hi):
    cf = os.path.join(CACHE, f'{pheno}_{gene}.txt')
    if os.path.exists(cf):
        rows = [l.split('\t') for l in open(cf).read().splitlines() if l]
    else:
        rows = []
        for _ in range(12):
            try:
                rows = [l.split('\t') for l in pysam.TabixFile(FF % pheno).fetch(str(chrom), int(lo), int(hi))]
                break
            except Exception:
                time.sleep(2)
        open(cf, 'w').write('\n'.join('\t'.join(r) for r in rows))
    return pd.DataFrame([{'snp': f[4].split(',')[0], 'ref': f[2], 'alt': f[3],
                          'beta': float(f[8]), 'se': float(f[9]), 'eaf': float(f[10])} for f in rows])

mr_rows = []
for g, (dsid, N, src) in PQTL.items():
    cr = COORD[g]; c19, p19 = cr['hg19'][0], cr['hg19'][1]; c38, p38 = cr['hg38'][0], cr['hg38'][1]
    p = og_cis(dsid, c19, p19 - WIN, p19 + WIN)
    if len(p) == 0:
        print(f'{g}: pQTL fetch EMPTY'); continue
    p = p[['rsid','ea','nea','beta','se','eaf','p','n']].rename(columns={'rsid':'snp'}).dropna().drop_duplicates('snp')
    gws = p[p['p'] < 5e-8]
    lead = gws.loc[gws['p'].idxmin()] if len(gws) else None
    n_gws = len(gws)
    for organ, (lab, pheno) in DISEASE.items():
        f = fg_cis(pheno, g, c38, p38 - WIN, p38 + WIN)
        m = p.merge(f, on='snp', suffixes=('.p', '.f'))
        m = m[m['se.f'] > 0].copy()
        m['by'] = np.where((m['alt']==m['ea'])&(m['ref']==m['nea']), m['beta.f'],
                   np.where((m['alt']==m['nea'])&(m['ref']==m['ea']), -m['beta.f'], np.nan))
        m = m[m['by'].notna()]
        if len(m):
            m.assign(gene=g, organ=organ).to_csv(f'{OUT}/{g}__{organ}.csv', index=False)
        if lead is not None:
            r = m[m['snp'] == lead['snp']]
            if len(r):
                bx, sx = lead['beta'], lead['se']; by = r['by'].iloc[0]; sy = r['se.f'].iloc[0]
                beta = by / bx; se = abs(sy / bx); pv = 2 * norm.sf(abs(beta / se))
                mr_rows.append({'gene':g,'source':src,'organ':organ,'disease':lab,'lead':lead['snp'],
                                'lead_p':lead['p'],'beta':beta,'OR':np.exp(beta),'p':pv})
    print(f'{g}: reg SNPs={len(p)}, GWS={n_gws}, lead={"-".join(lead[["snp"]].tolist()) if lead is not None else "none"}')

R = pd.DataFrame(mr_rows).sort_values(['gene','p'])
R.to_csv('out/screen_mr_results.csv', index=False)
print('\n=== lead-SNP Wald MR ===')
print(R.to_string(index=False, float_format=lambda x: '%.3g' % x))
