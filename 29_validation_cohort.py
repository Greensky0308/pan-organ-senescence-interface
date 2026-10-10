import gzip
import pandas as pd

# ---- GPL10687 probe->symbol ----
sym = {}
with gzip.open('cache/GPL10687.soft.gz', 'rt', errors='ignore') as f:
    for l in f:
        if l.startswith('!platform_table_begin'):
            next(f)
            for row in f:
                if row.startswith('!platform_table_end'): break
                p = row.rstrip('\n').split('\t')
                if len(p) > 2 and p[2].strip():
                    sym[p[0].strip()] = p[2].split('///')[0].strip()
            break
print('probe->symbol:', len(sym))

# ---- GSE25097 series matrix ----
tbl, meta, in_t = [], {}, False
with gzip.open('cache/GSE25097.gz', 'rt', errors='ignore') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'): in_t = True; continue
        if line.startswith('!series_matrix_table_end'): break
        if line.startswith('!Sample_'):
            k = line.split('\t', 1)[0].replace('!Sample_', '')
            meta.setdefault(k, []).append([v.strip('"\n') for v in line.rstrip('\n').split('\t')[1:]])
        if in_t: tbl.append(line.rstrip('\n').split('\t'))
hdr = [c.strip('"') for c in tbl[0][1:]]
mat = pd.DataFrame([[x.strip('"') for x in r[1:]] for r in tbl[1:]]).apply(pd.to_numeric, errors='coerce')
mat.index = [r[0].strip('"') for r in tbl[1:]]; mat.columns = hdr

tissue = None
for k, v in meta.items():
    for r in v:
        if r and any('cirrhotic liver' in x or 'healthy liver' in x for x in r): tissue = r; break
    if tissue: break
gm = dict(zip(hdr, tissue))
def grp(v):
    if 'cirrhotic liver' in v: return 'case'
    if 'healthy liver' in v: return 'control'
    return None
cols = [c for c in mat.columns if grp(gm.get(c)) is not None]
mat = mat[cols]
keep = [i for i in mat.index if i in sym]
e = mat.loc[keep].copy(); e.index = [sym[i] for i in keep]; e = e.groupby(level=0).mean()
e.to_csv('out/valid/liver_cirrhosis_expr.tsv', sep='\t')
pd.DataFrame({'sample': cols, 'group': [grp(gm[c]) for c in cols],
              'stratum': [grp(gm[c]) for c in cols]}).to_csv('out/valid/liver_cirrhosis_groups.tsv', sep='\t', index=False)
print('GSE25097 liver_cirrhosis:', e.shape, pd.Series([grp(gm[c]) for c in cols]).value_counts().to_dict())
