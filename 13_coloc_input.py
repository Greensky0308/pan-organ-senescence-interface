import os, json, time
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import pandas as pd, requests, pysam

TOK = os.environ.get('OPENGWAS_JWT', '')
H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}
OUT = 'out/coloc_input'; os.makedirs(OUT, exist_ok=True)
PQTL = {'MMP7': 'ebi-a-GCST90012056', 'MMP12': 'ebi-a-GCST90012070'}
DISEASE = {'kidney': ('DKD','DM_NEPHROPATHY'), 'liver': ('MASLD','NAFLD'),
           'heart': ('DCM','I9_NONISCHCARDMYOP'), 'lung': ('IPF','IPF'),
           'brain': ('AD','G6_ALZHEIMER'), 'colon': ('IBD','K11_IBD_STRICT'),
           'skin': ('Psoriasis','L12_PSORIASIS'), 'artery': ('Athero','I9_ATHSCLE'),
           'esophagus': ('EAC','C3_OESOPHAGUS_EXALLC')}
FF = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz'

phenos = {x['phenocode']: x for x in json.load(open('cache/ff_phenos.json'))}

# ---- pQTL cis (b37 11:102.3-103.0Mb) ----
s = requests.Session(); s.proxies = ({'https': os.environ['HTTPS_PROXY']} if os.environ.get('HTTPS_PROXY') else None)
for g, ds in PQTL.items():
    for _ in range(6):
        try:
            j = s.post('https://api.opengwas.io/api/associations', headers=H,
                       data=json.dumps({'id': [ds], 'variant': ['11:102300000-103000000']}), timeout=120).json()
            break
        except Exception:
            time.sleep(2)
    df = pd.DataFrame(j).rename(columns={'rsid': 'snp', 'beta': 'beta', 'se': 'se', 'eaf': 'eaf'})
    df = df[['snp', 'beta', 'se', 'eaf', 'n']].dropna()
    df['MAF'] = df['eaf'].clip(upper=0.5)
    df.to_csv(f'{OUT}/pqtl_{g}.csv', index=False)
    print(f'pQTL {g}: {len(df)} SNPs, N={df["n"].iloc[0]}')

# ---- FinnGen cis (hg38) ----
def fetch(url, tries=10):
    for _ in range(tries):
        try:
            return list(pysam.TabixFile(url).fetch('11', 102_300_000, 103_000_000))
        except Exception:
            time.sleep(2)
    return []

for organ, (lab, pheno) in DISEASE.items():
    rows = []
    for line in fetch(FF % pheno):
        f = line.split('\t')
        rows.append({'snp': f[4].split(',')[0], 'ref': f[2], 'alt': f[3],
                     'position': int(f[1]), 'beta': float(f[8]), 'se': float(f[9]),
                     'eaf': float(f[10])})
    df = pd.DataFrame(rows).dropna()
    df['MAF'] = df['eaf'].clip(upper=0.5)
    p = phenos[pheno]; N = p['num_cases'] + p['num_controls']; sc = p['num_cases'] / N
    df['N'] = N; df['s'] = sc
    df.to_csv(f'{OUT}/fg_{organ}.csv', index=False)
    print(f'FinnGen {organ:9s} {lab:9s}: {len(df)} SNPs, N={N}, s={sc:.3f}')
