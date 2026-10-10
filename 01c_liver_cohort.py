import gzip
import pandas as pd
RAW = 'raw'

def read_sm(path):
    tbl, meta, in_t = [], {}, False
    for line in gzip.open(path, 'rt'):
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
    for line in gzip.open(path, 'rt'):
        if line[0] in '^!#': continue
        p = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = p; ci = hdr.index('ID'); si = next(i for i, h in enumerate(hdr) if 'Gene symbol' in h); continue
        if len(p) <= max(ci, si): continue
        sym = p[si].split('///')[0].strip()
        if sym and sym not in ('---', ''): m[p[ci]] = sym
    return m

mat, hdr, meta = read_sm(f'{RAW}/GSE6764.sm.gz')
annot = read_annot(f'{RAW}/GPL570.annot.gz')
keep = [i for i in mat.index if i in annot]
e = mat.loc[keep].copy(); e.index = [annot[i] for i in keep]; e = e.groupby(level=0).mean()

raw = None
for k, v in meta.items():
    for r in v:
        if r and 'liver tissue' in r[0].lower(): raw = r; break
    if raw: break
gm = dict(zip(hdr, raw))
def stra(v):
    if 'cirrhotic' in v: return 'cirrhosis'
    if 'normal liver tissue' in v: return 'control'
    return None
cols = [c for c in e.columns if stra(gm.get(c)) is not None]
e = e[cols]
grp = pd.DataFrame({'sample': cols, 'group': ['case' if stra(gm[c]) == 'cirrhosis' else 'control' for c in cols],
                    'stratum': [stra(gm[c]) for c in cols]})
e.to_csv('out/cohorts/liver_expr.tsv', sep='\t')
grp.to_csv('out/cohorts/liver_groups.tsv', sep='\t', index=False)
print('liver(GSE6764 cirrhosis):', e.shape, grp['stratum'].value_counts().to_dict())
