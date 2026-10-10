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
                hdr = p; ci = hdr.index('ID')
                si = next(i for i, h in enumerate(hdr) if 'Gene symbol' in h); continue
            if len(p) <= max(ci, si): continue
            sym = p[si].split('///')[0].strip()
            if sym and sym not in ('---', ''): m[p[ci]] = sym
    return m

def row_with(meta, needle):
    for k, v in meta.items():
        for r in v:
            if r and any(needle.lower() in x.lower() for x in r): return r
    return None

def build(organ, smf, annotf, needle, ctrl, case, extra=None):
    mat, hdr, meta = read_sm(f'{RAW}/{smf}')
    annot = read_annot(f'{RAW}/{annotf}')
    keep = [i for i in mat.index if i in annot]
    e = mat.loc[keep].copy(); e.index = [annot[i] for i in keep]; e = e.groupby(level=0).mean()
    row = row_with(meta, needle)
    gm = dict(zip(hdr, [v.split(': ')[-1] for v in row]))
    ex = dict(zip(hdr, [v.split(': ')[-1] for v in row_with(meta, extra)])) if extra else None
    cols = [c for c in e.columns if case(gm.get(c)) is not None and (extra is None or ex.get(c) == 'Colon')]
    e = e[cols]
    g = [case(gm[c]) for c in cols]
    e.to_csv(f'{OUT}/{organ}_expr.tsv', sep='\t')
    pd.DataFrame({'sample': cols, 'group': g, 'stratum': g}).to_csv(f'{OUT}/{organ}_groups.tsv', sep='\t', index=False)
    print(f'{organ:9s} {e.shape}  ' + str(pd.Series(g).value_counts().to_dict()))

build('heart', 'GSE79962.sm.gz', 'GPL6244.annot.gz', 'Cardiac tissue',
      None, lambda v: 'case' if 'DCM' in v else ('control' if 'donor heart' in v else None))
# lung: IPF vs healthy
build('lung', 'GSE24206.sm.gz', 'GPL570.annot.gz', 'lung tissue',
      None, lambda v: 'case' if 'IPF' in v else ('control' if 'healthy' in v else None))
build('colon', 'GSE16879.sm.gz', 'GPL570.annot.gz', 'disease',
      None, lambda v: 'control' if v == 'Control' else ('case' if v in ('UC', 'CD') else None), extra='tissue')
build('skin', 'GSE14905.sm.gz', 'GPL570.annot.gz', 'Skin biopsy',
      None, lambda v: 'case' if ('lesional' in v and 'non-' not in v) else ('control' if v.endswith('normal') else None))
build('esophagus', 'GSE26886.sm.gz', 'GPL570.annot.gz', 'esophageal',
      None, lambda v: 'case' if 'adenocarcinoma' in v else ('control' if 'squamous epithelium' in v else None))
build('artery', 'GSE28829.sm.gz', 'GPL570.annot.gz', 'plaque',
      None, lambda v: 'case' if 'advanced' in v else ('control' if 'early' in v else None))

build('brain', 'GSE28146.sm.gz', 'GPL570.annot.gz', 'disease status',
      None, lambda v: 'control' if v == 'Control' else ('case' if v in ('Incipient', 'Moderate', 'Severe') else None))
