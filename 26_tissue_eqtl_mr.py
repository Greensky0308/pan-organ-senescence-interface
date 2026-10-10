import os, time, re
_pex = os.environ.get('HTTPS_PROXY')
if _pex:
    os.environ['https_proxy'] = os.environ['http_proxy'] = _pex
import numpy as np, pandas as pd, pysam
from scipy.stats import norm

ARC = 'out/causal'
BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz'
TM = {'kidney': ['Kidney_Cortex'], 'liver': ['Liver'], 'heart': ['Heart_Left_Ventricle', 'Heart_Atrial_Appendage'],
      'lung': ['Lung'], 'colon': ['Colon_Sigmoid', 'Colon_Transverse'],
      'skin': ['Skin_Sun_Exposed_Lower_leg', 'Skin_Not_Sun_Exposed_Suprapubic'],
      'artery': ['Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial'],
      'esophagus': ['Esophagus_Mucosa', 'Esophagus_Muscularis'],
      'brain': ['Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus']}
DIS = {'kidney': ('DKD', 'DM_NEPHROPATHY'), 'liver': ('Cirrhosis', 'CIRRHOSIS_BROAD'),
       'heart': ('DCM', 'I9_NONISCHCARDMYOP'), 'lung': ('IPF', 'IPF'), 'brain': ('AD', 'G6_ALZHEIMER'),
       'colon': ('IBD', 'K11_IBD_STRICT'), 'skin': ('Psoriasis', 'L12_PSORIASIS'),
       'artery': ('Athero', 'I9_ATHSCLE'), 'esophagus': ('EAC', 'C3_OESOPHAGUS_EXALLC')}

E = pd.read_csv(f'{ARC}/gtex_core_tissue_eqtl.csv')
rows = []
for organ, tissues in TM.items():
    lab, pheno = DIS[organ]
    task = []
    for g in E['gene'].unique():
        sub = E[(E['gene'] == g) & (E['tissue'].isin(tissues))]
        if sub.empty:
            continue
        lead = sub.loc[sub['pval'].idxmin()]
        m = re.match(r'(chr[\dXY]+)_(\d+)_([A-Za-z]+)_([A-Za-z]+)_b38', lead['variant_id'])
        if not m:
            continue
        task.append((g, m.group(1).replace('chr', ''), int(m.group(2)), m.group(3), m.group(4),
                     lead['variant_id'], lead['tissue'], lead['pval'], lead['slope']))
    if not task:
        continue
    for chrom in sorted(set(t[1] for t in task)):
        grp = [t for t in task if t[1] == chrom]
        lo, hi = min(t[2] for t in grp), max(t[2] for t in grp)
        lines = []
        for _ in range(6):
            try:
                lines = list(pysam.TabixFile(BASE % pheno).fetch(chrom, lo - 1, hi)); break
            except Exception:
                time.sleep(2)
        idx = {int(f.split('\t')[1]): f.split('\t') for f in lines}
        for (g, c, pos, ref, alt, vid, tis, ep, slope) in grp:
            f = idx.get(pos)
            if not f:
                continue
            fref, falt, by, sy = f[2], f[3], float(f[8]), float(f[9])
            if falt == alt and fref == ref:
                b = by
            elif falt == ref and fref == alt:
                b = -by
            else:
                continue
            beta = b / slope; se = abs(sy / slope)
            rows.append({'gene': g, 'organ': organ, 'disease': lab, 'tissue': tis, 'snp': vid,
                         'eqtl_p': ep, 'beta': beta, 'OR': float(np.exp(beta)),
                         'p': float(2 * norm.sf(abs(beta / se)))})
R = pd.DataFrame(rows).sort_values('p')
R.to_csv(f'{ARC}/eqtl_tissue_mr.csv', index=False)
print(R.to_string(index=False, float_format=lambda x: '%.4g' % x))
print('rows:', len(R))
