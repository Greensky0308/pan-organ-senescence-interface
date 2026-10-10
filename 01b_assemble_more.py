import gzip, os
import numpy as np, pandas as pd

RAW = 'raw'; OUT = 'out/cohorts'

def read_sm(path):
    tbl, in_t, meta = [], False, {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!series_matrix_table_begin'): in_t = True; continue
            if line.startswith('!series_matrix_table_end'): break
            if line.startswith('!Sample_'):
                k = line.split('\t', 1)[0].replace('!Sample_', '')
                meta.setdefault(k, []).append([v.strip('"\n') for v in line.rstrip('\n').split('\t')[1:]])
            if in_t: tbl.append(line.rstrip('\n').split('\t'))
    if len(tbl) > 1:
        hdr = [c.strip('"') for c in tbl[0][1:]]
        mat = pd.DataFrame([[x.strip('"') for x in r[1:]] for r in tbl[1:]]).apply(pd.to_numeric, errors='coerce')
        mat.index = [r[0].strip('"') for r in tbl[1:]]; mat.columns = hdr
    else:
        hdr = meta.get('geo_accession', [[]])[0]; mat = pd.DataFrame(columns=hdr)
    md = pd.DataFrame({k: v[0] for k, v in meta.items() if len(v) == 1}, index=hdr)
    return mat, md

def read_annot(path):
    m = {}; hdr = None
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] in '^!#' or line.startswith('#!'): continue
            p = line.rstrip('\n').split('\t')
            if hdr is None:
                hdr = p; ci = hdr.index('ID')
                si = next((i for i, h in enumerate(hdr) if 'Gene symbol' in h or h == 'Symbol'), None)
                continue
            if len(p) > max(ci, si):
                s = p[si].split('///')[0].strip()
                if s and s != '---': m[p[ci]] = s
    return m

def collapse(mat, id2sym):
    keep = [i for i in mat.index if i in id2sym]
    m = mat.loc[keep].copy(); m.index = [id2sym[i] for i in keep]
    return m.groupby(level=0).mean()

def char_row(path, field):
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('!Sample_characteristics_ch1'):
                parts = [p.strip('"') for p in line.rstrip('\n').split('\t')[1:]]
                ne = next((p for p in parts if p), '')
                if field in ne:
                    return [(p.split(':', 1)[1].strip() if ':' in p else p) for p in parts]
    return None

def save(organ, expr, groups):
    expr.to_csv(f'{OUT}/{organ}_expr.tsv', sep='\t')
    pd.DataFrame({'sample': expr.columns, 'group': groups,
                  'stratum': groups}).to_csv(f'{OUT}/{organ}_groups.tsv', sep='\t', index=False)
    print(f'{organ:11} {expr.shape} {pd.Series(groups).value_counts().to_dict()}')

import re
mat, md = read_sm(f'{RAW}/GSE5281.sm.gz')
expr = collapse(mat, read_annot(f'{RAW}/GPL570.annot.gz'))
rows = []
with gzip.open(f'{RAW}/GSE5281.sm.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('!Sample_characteristics_ch1'):
            rows.append([p.strip('"') for p in line.rstrip('\n').split('\t')[1:]])
dx = []
for i in range(len(expr.columns)):
    v = ''
    for r in rows:
        if i < len(r) and re.search(r'disease state', r[i], re.I):
            v = r[i].split(':', 1)[1].strip(); break
    dx.append('case' if 'alzheimer' in v.lower() else ('control' if 'normal' in v.lower() else 'other'))
save('brain', expr, dx)

mat, md = read_sm(f'{RAW}/GSE43292.sm.gz')
expr = collapse(mat, read_annot(f'{RAW}/GPL6244.annot.gz'))
v = char_row(f'{RAW}/GSE43292.sm.gz', 'tissue')
save('artery', expr, ['case' if 'Atheroma' in x else 'control' for x in v])

mat, md = read_sm(f'{RAW}/GSE13898.sm.gz')
expr = collapse(mat, read_annot(f'{RAW}/GPL6102.annot.gz'))
v = char_row(f'{RAW}/GSE13898.sm.gz', 'Pathology')
save('esophagus', expr, ['case' if x == 'adenocarcinoma' else ('control' if x == 'normal' else 'other') for x in v])

mat, md = read_sm(f'{RAW}/GSE20950.sm.gz')
expr = collapse(mat, read_annot(f'{RAW}/GPL570.annot.gz'))
v = char_row(f'{RAW}/GSE20950.sm.gz', 'metabolic status')
save('adipose', expr, ['case' if 'resistant' in x else 'control' for x in v])

mat, md = read_sm(f'{RAW}/GSE7305.sm.gz')
expr = collapse(mat, read_annot(f'{RAW}/GPL570.annot.gz'))
t = md['title'].tolist()
save('endometrium', expr, ['control' if 'Normal' in x else 'case' for x in t])
print('DONE')
