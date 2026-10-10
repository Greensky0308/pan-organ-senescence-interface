# Build the cohort table and the cross-organ core-gene table (.xlsx + .docx).
import os
import numpy as np, pandas as pd
from docx import Document
from docx.shared import Pt, Inches

R = 'out'; MAIN = 'out/figures'; SUPP = 'out/tables'
os.makedirs(MAIN, exist_ok=True); os.makedirs(SUPP, exist_ok=True)

ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
LBL = {'kidney': 'Kidney', 'liver': 'Liver', 'heart': 'Heart', 'lung': 'Lung', 'colon': 'Colon',
       'skin': 'Skin', 'artery': 'Artery', 'esophagus': 'Esophagus', 'brain': 'Brain'}
DISEASE = {'kidney': 'Diabetic kidney disease', 'liver': 'Cirrhosis', 'heart': 'Dilated cardiomyopathy',
           'lung': 'Idiopathic pulmonary fibrosis', 'colon': 'Inflammatory bowel disease', 'skin': 'Psoriasis',
           'artery': 'Atherosclerosis', 'esophagus': 'Esophageal adenocarcinoma', 'brain': 'Alzheimer disease'}
COHORT = {'kidney': 'GSE104954', 'liver': 'GSE6764', 'heart': 'GSE141910', 'lung': 'GSE150910', 'colon': 'GSE75214',
          'skin': 'GSE13355', 'artery': 'GSE43292', 'esophagus': 'GSE13898', 'brain': 'GSE5281'}
PLAT = {'kidney': 'Microarray', 'liver': 'Microarray', 'heart': 'RNA-seq', 'lung': 'RNA-seq', 'colon': 'Microarray',
        'skin': 'Microarray', 'artery': 'Microarray', 'esophagus': 'Microarray', 'brain': 'Microarray'}

en = pd.read_csv(f'{R}/interface_enrichment.csv').set_index('organ').loc[ORGS]
M = pd.read_csv(f'{R}/interface_matrix.csv', index_col=0)
AGE = pd.read_csv(f'{R}/gtex_age_summary.csv').set_index('gene')
TAU = pd.read_csv(f'{R}/interface_gtex_tau.csv', index_col=0)['tau']
ms = pd.read_csv(f'{R}/multisig_enrichment.csv')
SIGORD = ['SenCat', 'SenMayo', 'Fridman', 'ReactomeSASP', 'Reactome', 'GO0090398']


DEF = {'kidney': ({'LD', 'TN'}, {'DN', 'FSGS', 'MCD', 'RPGN'}), 'liver': ({'control'}, {'cirrhosis'}),
       'heart': ({'Non-Failing Donor'}, {'Dilated cardiomyopathy (DCM)'}), 'lung': ({'control'}, {'ipf'}),
       'colon': ({'control'}, {'ulcerative colitis', 'CD'}), 'skin': ({'control'}, {'case'}),
       'artery': ({'control'}, {'case'}), 'esophagus': ({'control'}, {'case'}), 'brain': ({'control'}, {'case'})}


def n_samples(o):
    g = pd.read_csv(f'{R}/cohorts/{o}_groups.tsv', sep='\t')
    ctrl, case = DEF[o]
    return int(g['stratum'].isin(case).sum()), int(g['stratum'].isin(ctrl).sum())


def write_table(df, path_xlsx, path_docx, caption=''):
    df.to_excel(path_xlsx, index=False)
    doc = Document()
    st = doc.styles['Normal']; st.font.name = 'Times New Roman'; st.font.size = Pt(9)
    if caption:
        p = doc.add_paragraph(); r = p.add_run(caption); r.bold = True; r.font.size = Pt(10); r.font.name = 'Times New Roman'
    t = doc.add_table(rows=1, cols=df.shape[1]); t.style = 'Table Grid'
    for j, c in enumerate(df.columns):
        cell = t.rows[0].cells[j]; cell.text = str(c)
        for pr in cell.paragraphs[0].runs: pr.bold = True; pr.font.size = Pt(9); pr.font.name = 'Times New Roman'
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, c in enumerate(df.columns):
            cells[j].text = str(row[c])
            for pr in cells[j].paragraphs[0].runs: pr.font.size = Pt(9); pr.font.name = 'Times New Roman'
    doc.save(path_docx)
    print('saved', path_xlsx, '+', os.path.basename(path_docx))


# ================= cohort overview =================
rows = []
for o in ORGS:
    nc, nk = n_samples(o)
    rows.append({'Organ': LBL[o], 'Disease': DISEASE[o], 'Cohort': COHORT[o], 'Platform': PLAT[o],
                 'Cases': nc, 'Controls': nk, 'Signature genes (up)': int(en.loc[o, 'n_sig']),
                 'SenCat overlap': int(en.loc[o, 'overlap']), 'Enrichment fold': round(float(en.loc[o, 'fold']), 2),
                 'P value': ('%.1e' % en.loc[o, 'p'])})
T1 = pd.DataFrame(rows)
write_table(T1, f'{MAIN}/Table1.xlsx', f'{MAIN}/Table1.docx',
            'Nine organ disease transcriptome cohorts and senescence-signature enrichment.')

# ================= shared core genes =================
sub = M[M['n_organs'] >= 4].copy().sort_values('n_organs', ascending=False)
rows = []
for g in sub.index:
    orgs = [LBL[o] for o in ORGS if sub.loc[g, o] == 1]
    rows.append({'Gene': g, 'Shared organs (n)': int(sub.loc[g, 'n_organs']), 'Organs': ', '.join(orgs),
                 'Age-up tissues (of 49)': (int(AGE.loc[g, 'n_sig_pos']) if g in AGE.index else ''),
                 'Mean age rho': (round(float(AGE.loc[g, 'mean_rho']), 3) if g in AGE.index else ''),
                 'GTEx tau': (round(float(TAU[g]), 3) if g in TAU.index and pd.notna(TAU.get(g)) else '')})
T2 = pd.DataFrame(rows)
write_table(T2, f'{MAIN}/Table2.xlsx', f'{MAIN}/Table2.docx',
            'Interface genes shared by >=4 organs (universal core).')
print('T1', T1.shape, 'T2', T2.shape)
