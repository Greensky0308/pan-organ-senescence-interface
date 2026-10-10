import os, json, time
PROXY = os.environ.get('HTTPS_PROXY')
if PROXY: os.environ['https_proxy'] = os.environ['http_proxy'] = PROXY
import pandas as pd, requests

TOK = os.environ.get('OPENGWAS_JWT', ''); H = {'Authorization': 'Bearer ' + TOK, 'Content-Type': 'application/json'}
S = requests.Session(); S.proxies = ({'https': os.environ['HTTPS_PROXY']} if os.environ.get('HTTPS_PROXY') else None)
OUT = 'out/coloc_input_lyz'; os.makedirs(OUT, exist_ok=True)
REG = '12:69400000-70100000'          # LYZ hg19 ±~300kb
VAL = {'ieu-a-294': 'IBD_Liu2015', 'ieu-a-292': 'IBD_Jostins2012', 'ieu-a-12': 'Crohn_Liu2015'}

def og(dsid, var):
    for _ in range(6):
        try:
            j = S.post('https://api.opengwas.io/api/associations', headers=H,
                       data=json.dumps({'id': [dsid], 'variant': [var]}), timeout=120).json()
            if isinstance(j, list) and j and 'beta' in j[0]:
                return pd.DataFrame(j)
        except Exception:
            time.sleep(2)
    return pd.DataFrame()

p = og('prot-a-1821', REG)[['rsid','ea','nea','beta','se','eaf','p','n']].dropna().drop_duplicates('rsid')
p.to_csv(f'{OUT}/pqtl_LYZ.csv', index=False)
print('pQTL LYZ snp:', len(p))
for dsid, lab in VAL.items():
    info = S.get(f'https://api.opengwas.io/api/gwasinfo/{dsid}', headers=H, timeout=60)
    g = og(dsid, REG)
    if len(g) == 0:
        print(lab, 'no data'); continue
    g = g[['rsid','ea','nea','beta','se','eaf','p']].dropna().drop_duplicates('rsid')
    try:
        meta = S.get(f'https://api.opengwas.io/api/gwasinfo', headers=H, params={'id': dsid}, timeout=60).json()
    except Exception:
        meta = {}
    g.to_csv(f'{OUT}/{lab}.csv', index=False)
    print(f'{lab}: gwas snps={len(g)}')
