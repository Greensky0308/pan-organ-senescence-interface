import gzip, glob, os, re
from collections import defaultdict
import numpy as np, pandas as pd

RAW = 'raw'
OUT = 'out/cohorts'
os.makedirs(OUT, exist_ok=True)
KID = os.environ.get('KID_DIR', 'data/kidney_prior')

def read_sm(path):
    rows, tbl, in_t, meta = [], [], False, defaultdict(list)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!series_matrix_table_begin'): in_t = True; continue
            if line.startswith('!series_matrix_table_end'): break
            if line.startswith('!Sample_'):
                k = line.split('\t', 1)[0].replace('!Sample_', '')
                vals = [v.strip('"\n') for v in line.rstrip('\n').split('\t')[1:]]
                meta[k].append(vals)
            if in_t: tbl.append(line.rstrip('\n').split('\t'))
    if len(tbl) > 1:
        hdr = [c.strip('"') for c in tbl[0][1:]]
        mat = pd.DataFrame([[x.strip('"') for x in r[1:]] for r in tbl[1:]]).apply(pd.to_numeric, errors='coerce')
        mat.index = [r[0].strip('"') for r in tbl[1:]]; mat.columns = hdr
    else:
        hdr = meta.get('geo_accession', [[]])[0]
        mat = pd.DataFrame(columns=hdr)
    md = pd.DataFrame({k: v[0] for k, v in meta.items() if len(v) == 1}, index=hdr)
    return mat, md

def read_annot(path):
    m = {}
    with gzip.open(path, 'rt') as fh:
        hdr = None
        for line in fh:
            if line.startswith('^') or line.startswith('!'): continue
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if hdr is None:
                hdr = p;
                ci = hdr.index('ID'); si = next((i for i, h in enumerate(hdr) if 'Gene symbol' in h or h == 'Symbol'), None)
                if si is None: raise ValueError('no symbol col: ' + str(hdr))
                continue
            if len(p) <= max(ci, si): continue
            sym = p[si].split('///')[0].strip()
            if sym and sym not in ('---', ''): m[p[ci]] = sym.split('///')[0].strip()
    return m

def char_row(path, field, ch='characteristics_ch1'):
    """Return the per-sample value of a characteristics row whose first value contains `field`. """
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!Sample_' + ch):
                parts = [p.strip('"') for p in line.rstrip('\n').split('\t')[1:]]
                if parts and field in parts[0]:
                    return [p.split(': ', 1)[-1] for p in parts]
    return None

_hg = pd.read_csv(f'{RAW}/hgnc.txt', sep='\t', dtype=str, low_memory=False)
ENS2SYM = dict(zip(_hg['ensembl_gene_id'].str.split('.').str[0], _hg['symbol']))

def ens_collapse(df):
    df = df.copy()
    df.index = [ENS2SYM.get(str(i).split('.')[0]) for i in df.index]
    df = df[df.index.notna()]
    return df.groupby(level=0).mean()

def collapse(mat, id2sym):
    keep = [i for i in mat.index if i in id2sym]
    m = mat.loc[keep].copy(); m.index = [id2sym[i] for i in keep]
    m = m.groupby(level=0).mean()
    return m

# ---------- brain (GSE33000, GPL4372) ----------
mat, md = read_sm(f'{RAW}/GSE33000.sm.gz')
sym = read_annot(f'{RAW}/GPL4372.annot.gz')
expr = collapse(mat, sym)
dxvals = char_row(f'{RAW}/GSE33000.sm.gz', 'disease status', 'characteristics_ch2')
grp = {s: ('AD' if 'Alzheimer' in v else ('HD' if 'Huntington' in v else 'control'))
       for s, v in zip(md.index, dxvals)}
expr.to_csv(f'{OUT}/brain_expr.tsv', sep='\t')
pd.DataFrame({'sample': list(grp), 'group': [grp[s] for s in grp], 'stratum': [grp[s] for s in grp]}).to_csv(f'{OUT}/brain_groups.tsv', sep='\t', index=False)
print('brain', expr.shape, pd.Series(grp).value_counts().to_dict())

# ---------- colon (GSE75214, GPL6244) ----------
mat, md = read_sm(f'{RAW}/GSE75214.sm.gz')
sym = read_annot(f'{RAW}/GPL6244.annot.gz')
expr = collapse(mat, sym)
g = pd.Series(char_row(f'{RAW}/GSE75214.sm.gz', 'disease'), index=expr.columns).replace({"Crohn's disease": 'CD'})
expr.to_csv(f'{OUT}/colon_expr.tsv', sep='\t')
pd.DataFrame({'sample': expr.columns, 'group': [('case' if x in ('ulcerative colitis', 'CD') else 'control') for x in g],
              'stratum': g.values}).to_csv(f'{OUT}/colon_groups.tsv', sep='\t', index=False)
print('colon', expr.shape, g.value_counts().to_dict())

# ---------- skin (GSE13355, GPL570) ----------
mat, md = read_sm(f'{RAW}/GSE13355.sm.gz')
sym = read_annot(f'{RAW}/GPL570.annot.gz')
expr = collapse(mat, sym)
vals = char_row(f'{RAW}/GSE13355.sm.gz', 'skin')
g = ['control' if 'normal' in v else ('case' if v.startswith('involved') else 'other') for v in vals]
expr.to_csv(f'{OUT}/skin_expr.tsv', sep='\t')
pd.DataFrame({'sample': expr.columns, 'group': g, 'stratum': g}).to_csv(f'{OUT}/skin_groups.tsv', sep='\t', index=False)
print('skin', expr.shape, pd.Series(g).value_counts().to_dict())

# ---------- lung (GSE150910, gene counts) ----------
cnt = pd.read_csv(f'{RAW}/GSE150910_counts.csv.gz', index_col=0)
cnt = cnt.groupby(level=0).sum()
mat, md = read_sm(f'{RAW}/GSE150910.sm.gz')
title2dx = {}
with gzip.open(f'{RAW}/GSE150910.sm.gz', 'rt') as fh:
    t = d = None
    for line in fh:
        if line.startswith('!Sample_title'):
            t = [p.strip('"') for p in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and 'diagnosis' in line:
            d = [p.strip('"').split(': ')[-1] for p in line.rstrip('\n').split('\t')[1:]]
    title2dx = dict(zip(t, d))
cols = [c for c in cnt.columns if title2dx.get(c) in ('ipf', 'control', 'chp')]
expr = cnt[cols]
g = [title2dx[c] for c in cols]
expr.to_csv(f'{OUT}/lung_expr.tsv', sep='\t')
pd.DataFrame({'sample': cols, 'group': [('case' if x in ('ipf', 'chp') else 'control') for x in g], 'stratum': g}).to_csv(f'{OUT}/lung_groups.tsv', sep='\t', index=False)
print('lung', expr.shape, pd.Series(g).value_counts().to_dict())

# ---------- liver (GSE135251, counts per-sample) ----------
mat, md = read_sm(f'{RAW}/GSE135251.sm.gz')
files = sorted(glob.glob(f'{RAW}/_p/GSM*_*.counts.txt.gz'))
if len(files) < 200:
    os.makedirs(f'{RAW}/_p', exist_ok=True)
    os.system(f'tar xf {RAW}/GSE135251_RAW.tar -C {RAW}/_p 2>/dev/null')
    files = sorted(glob.glob(f'{RAW}/_p/GSM*_*.counts.txt.gz'))
series = {}
for p in files:
    gsm = os.path.basename(p).split('_')[0]
    d = {}
    try:
        with gzip.open(p, 'rt') as fh:
            for l in fh:
                if '\t' not in l: continue
                a, b = l.rstrip('\n').split('\t')[:2]
                try: d[a] = float(b)
                except ValueError: pass
    except EOFError:
        print('  skip truncated:', os.path.basename(p)); continue
    series[gsm] = d
liv = ens_collapse(pd.DataFrame(series))
grp_paper = dict(zip(md.index, char_row(f'{RAW}/GSE135251.sm.gz', 'group in paper')))
liv = liv[[c for c in liv.columns if c in grp_paper]]
liv.to_csv(f'{OUT}/liver_expr.tsv', sep='\t')
st = [grp_paper[c] for c in liv.columns]
pd.DataFrame({'sample': liv.columns, 'group': [('control' if x == 'control' else 'case') for x in st],
              'stratum': st}).to_csv(f'{OUT}/liver_groups.tsv', sep='\t', index=False)
print('liver', liv.shape, pd.Series(st).value_counts().to_dict())

# ---------- heart (GSE141910, normalized per-sample) ----------
mat, md = read_sm(f'{RAW}/GSE141910.sm.gz')
files = sorted(glob.glob(f'{RAW}/_h/GSM*.csv.gz'))
if len(files) < 300:
    os.makedirs(f'{RAW}/_h', exist_ok=True)
    os.system(f'tar xf {RAW}/GSE141910_RAW.tar -C {RAW}/_h 2>/dev/null')
    files = sorted(glob.glob(f'{RAW}/_h/GSM*.csv.gz'))
series = {}
for p in files:
    gsm = os.path.basename(p).split('_')[0]
    df = pd.read_csv(p, index_col=0).iloc[:, 0]
    series[gsm] = df
hea = ens_collapse(pd.DataFrame(series))
etio = dict(zip(md.index, char_row(f'{RAW}/GSE141910.sm.gz', 'etiology')))
hea = hea[[c for c in hea.columns if c in etio]]
hea.to_csv(f'{OUT}/heart_expr.tsv', sep='\t')
st = ['control' if 'Non-Failing' in etio[c] else 'case' for c in hea.columns]
pd.DataFrame({'sample': hea.columns, 'group': st, 'stratum': [etio[c] for c in hea.columns]}).to_csv(f'{OUT}/heart_groups.tsv', sep='\t', index=False)
print('heart', hea.shape, pd.Series([etio[c] for c in hea.columns]).value_counts().to_dict())

gi = pd.read_csv(f'{KID}/gene_info.gz', sep='\t', dtype=str)
gi = gi[gi['#tax_id'] == '9606'][['GeneID', 'Symbol']].drop_duplicates('GeneID')
eg2sym = dict(zip(gi['GeneID'], gi['Symbol']))
k = pd.read_csv(f'{KID}/expr.tsv.gz', sep='\t', index_col=0)
k.index = [eg2sym.get(i.rsplit('_', 1)[0], None) for i in k.index]
k = k[k.index.notna()].groupby(level=0).mean()
grp = pd.read_csv(f'{KID}/grp.txt', header=None)[0].tolist()
k.to_csv(f'{OUT}/kidney_expr.tsv', sep='\t')
pd.DataFrame({'sample': k.columns, 'group': [('control' if g in ('LD', 'TN') else 'case') for g in grp],
              'stratum': grp}).to_csv(f'{OUT}/kidney_groups.tsv', sep='\t', index=False)
print('kidney', k.shape, pd.Series(grp).value_counts().to_dict())
print('\nDONE')
