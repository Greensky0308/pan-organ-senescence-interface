import json, subprocess, tempfile, os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
M = pd.read_csv('out/interface_matrix.csv', index_col=0)
SHARED = list(M[M['n_organs'] >= 3].index)

print('=== Interface module score (case vs control) ===')
res = []
for o in ORGS:
    expr = pd.read_csv(f'out/cohorts/{o}_expr.tsv', sep='\t', index_col=0)
    grp = pd.read_csv(f'out/cohorts/{o}_groups.tsv', sep='\t')
    gmap = dict(zip(grp['sample'], grp['group']))
    genes = [g for g in M.index if M.loc[g, o] == 1 and g in expr.index]
    sub = expr.loc[genes]
    z = sub.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) + 1e-9), axis=1)
    score = z.mean(axis=0)
    gr = pd.Series({s: gmap.get(s, 'na') for s in expr.columns})
    case = score[gr == 'case']; ctl = score[gr == 'control']
    p = mannwhitneyu(case, ctl, alternative='two-sided')[1]
    res.append({'organ': o, 'n_genes': len(genes), 'case_med': round(case.median(), 2),
                'ctrl_med': round(ctl.median(), 2), 'p': p})
    print('%-7s genes=%3d | case %.2f vs ctrl %.2f | p=%.3g' % (o, len(genes), case.median(), ctl.median(), p))
pd.DataFrame(res).to_csv('out/interface_score_per_organ.csv', index=False)

# ---------- 2. gprofiler ----------
def gprof(genes, name):
    if len(genes) < 5:
        print(name, 'too few'); return None
    body = {'organism': 'hsapiens', 'query': genes, 'sources': ['GO:BP', 'KEGG', 'REAC'],
            'user_threshold': 0.05, 'no_evidences': False}
    f = tempfile.NamedTemporaryFile('w', suffix='.json', delete=False); json.dump(body, f); f.close()
    out = subprocess.run(['curl', '-s', '-m', '120', '-x', os.environ.get('HTTPS_PROXY', ''),
                          '-X', 'POST', 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                          '-H', 'Content-Type: application/json', '--data', '@' + f.name],
                         capture_output=True, text=True).stdout
    os.unlink(f.name)
    try: r = json.loads(out)['result']
    except Exception: print(name, 'parse fail'); return None
    def gnames(x):
        out = []
        for it in (x.get('intersections') or []):
            if isinstance(it, (list, tuple)):
                if it: out.append(it[0])
            elif it:
                out.append(it)
        return ','.join(map(str, out))
    df = pd.DataFrame([{'term': x['name'], 'source': x['source'], 'p': x['p_value'],
                        'n': x.get('intersection_size'), 'genes': gnames(x)} for x in r])
    df.to_csv(f'out/pathway_{name}.csv', index=False)
    print(f'\n--- {name} ({len(genes)} genes): top terms ---')
    if len(df): print(df.sort_values('p').head(12).to_string(index=False)[:1400])
    return df

gprof(SHARED, 'shared_ge3organs')
gprof(list(M.index), 'all_interface')
gprof(list(M[M['n_organs'] == 1].index), 'organ_specific')

hg = pd.read_csv('raw/hgnc.txt', sep='\t', dtype=str, low_memory=False)
ENS2SYM = dict(zip(hg['ensembl_gene_id'].str.split('.').str[0], hg['symbol']))
g = pd.read_csv('raw/gtex_median.gct.gz', sep='\t', skiprows=2)
g['sym'] = [ENS2SYM.get(str(n).split('.')[0]) for n in g['Name']]
g = g[g['sym'].notna()].groupby('sym').mean(numeric_only=True)
tpm = g.values
mx = np.maximum(tpm.max(1, keepdims=True), 1e-9)
tau = (1 - tpm / mx).sum(1) / (tpm.shape[1] - 1)
T = pd.DataFrame({'n_tissues_tpm1': (tpm > 1).sum(1), 'tau': tau}, index=g.index).reindex(M.index)
T.to_csv('out/interface_gtex_tau.csv')
print('\n=== GTEx tau (closer to 1 = more tissue-specific) ===')
for lab, m in [('1 organ', M['n_organs'] == 1), ('2 organs', M['n_organs'] == 2), ('>=3 organs', M['n_organs'] >= 3)]:
    sub = T[m.values]
    print('%-9s n=%3d tau median=%.3f | median #tissues=%.0f' % (lab, sub['tau'].notna().sum(),
          sub['tau'].median(), sub['n_tissues_tpm1'].median()))
