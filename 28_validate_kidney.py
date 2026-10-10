import gzip, os
import pandas as pd

RAW = 'raw'; OUT = 'out/valid'; os.makedirs(OUT, exist_ok=True)

def read_sm(path):
    tbl, meta, in_t = [], {}, False
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!series_matrix_table_begin'): in_t = True; continue
            if line.startswith('!series_matrix_table_end'): break
            if line.startswith('!Sample_'):
                k = line.split('\t', 1)[0].replace('!Sample_', '')
                meta.setdefault(k, []).append([v.strip('"\n') for v in line.rstrip('\n').split('\t')[1:]])
            if in_t: tbl.append(line.rstrip('\n').split('\t'))
    hdr = [c.strip('"') for c in tbl[0][1:]]
    mat = pd.DataFrame([[x.strip('"') for x in r[1:]] for r in tbl[1:]]).apply(pd.to_numeric, errors='coerce')
    mat.index = [r[0].strip('"') for r in tbl[1:]]; mat.columns = hdr
    return mat, hdr, meta

def read_annot(path):
    m = {}; hdr = None
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] in '^!#': continue
            p = line.rstrip('\n').split('\t')
            if hdr is None:
                hdr = p; ci = hdr.index('ID'); si = next(i for i, h in enumerate(hdr) if 'Gene symbol' in h); continue
            if len(p) <= max(ci, si): continue
            s = p[si].split('///')[0].strip()
            if s and s not in ('---', ''): m[p[ci]] = s
    return m

def row_vals(meta, needle):
    for k, v in meta.items():
        for r in v:
            if r and any(needle.lower() in x.lower() for x in r):
                return [x.split(': ')[-1] for x in r]
    return None

g571 = {}
with open('cache/gpl571.soft') as fh:
    started = False
    for line in fh:
        if line.startswith('!platform_table_begin'): started = True; continue
        if not started: continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 11: continue
        s = p[10].split('///')[0].strip()
        if s and s not in ('---', '', 'Gene Symbol'): g571[p[0]] = s
print('GPL571 probes mapped:', len(g571))
mat, hdr, meta = read_sm(f'{RAW}/GSE30529.sm.gz')
keep = [i for i in mat.index if i in g571]
e = mat.loc[keep].copy(); e.index = [g571[i] for i in keep]; e = e.groupby(level=0).mean()
vals = None
for k, v in meta.items():
    for r in v:
        if r and any('Tubuli' in x for x in r): vals = [x.split(': ')[-1] for x in r]
gm = dict(zip(hdr, vals))
cols = [c for c in e.columns if 'Tubuli' in gm.get(c, '')]
e = e[cols]; g = ['control' if 'control' in gm[c] else 'case' for c in cols]
e.to_csv(f'{OUT}/kidney_expr.tsv', sep='\t')
pd.DataFrame({'sample': cols, 'group': g, 'stratum': g}).to_csv(f'{OUT}/kidney_groups.tsv', sep='\t', index=False)
print('kidney', e.shape, pd.Series(g).value_counts().to_dict())
