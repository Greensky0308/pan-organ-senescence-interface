import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

R = 'out'; OUT = 'out/figures'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'], 'font.size': 7,
                     'axes.linewidth': 0.6, 'axes.labelsize': 6.2, 'xtick.labelsize': 6.2, 'ytick.labelsize': 6.2,
                     'legend.fontsize': 6.2, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5})
CB = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#999999']
ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
LBL = {'kidney': 'DKD', 'liver': 'Cirrhosis', 'heart': 'DCM', 'lung': 'IPF', 'colon': 'IBD',
       'skin': 'Psoriasis', 'artery': 'Athero', 'esophagus': 'EAC', 'brain': 'AD'}
labs = [LBL[o] for o in ORGS]
M = pd.read_csv(f'{R}/interface_matrix.csv', index_col=0)
AGE = pd.read_csv(f'{R}/gtex_age_summary.csv')
en = pd.read_csv(f'{R}/interface_enrichment.csv').set_index('organ').loc[ORGS]


def despine(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
def plab(ax, L):
    ax.text(-0.16, 1.05, L, transform=ax.transAxes, fontsize=8, va='bottom', ha='left')
def save(fig, name):
    fig.savefig(f'{OUT}/{name}.pdf', dpi=600)
    fig.savefig(f'{OUT}/{name}.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig); print('saved', name)


fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.55, wspace=0.32, left=0.09, right=0.975, top=0.91, bottom=0.10)
gsa = gs[0, 0].subgridspec(1, 2, wspace=0.42)
A1 = fig.add_subplot(gsa[0]); A2 = fig.add_subplot(gsa[1])
B = fig.add_subplot(gs[0, 1]); C = fig.add_subplot(gs[1, 0]); D = fig.add_subplot(gs[1, 1])
for i, (lab, m) in enumerate([('background', AGE[~AGE['is_iface']]), ('interface', AGE[AGE['is_iface']])]):
    v = m['mean_rho'].values.astype(float); gx = np.linspace(v.min(), v.max(), 200)
    A1.plot(gx, gaussian_kde(v)(gx), color=CB[i], lw=1.2, label='%s (n=%d)' % (lab, len(m)))
    A1.fill_between(gx, gaussian_kde(v)(gx), color=CB[i], alpha=0.18)
A1.set_xlabel('mean age-corr rho'); A1.set_ylabel('density'); A1.set_ylim(0, A1.get_ylim()[1] * 1.32)
for i, (lab, m) in enumerate([('background', AGE[~AGE['is_iface']]), ('interface', AGE[AGE['is_iface']])]):
    v = m['n_sig_pos'].values.astype(float); gx = np.linspace(-2, 27, 200)
    A2.plot(gx, gaussian_kde(v, bw_method=0.35)(gx), color=CB[i], lw=1.2, label='%s (n=%d)' % (lab, len(m)))
    A2.fill_between(gx, gaussian_kde(v, bw_method=0.35)(gx), color=CB[i], alpha=0.18)
A2.set_xlabel('# tissues with age↑'); A2.set_ylim(0, A2.get_ylim()[1] * 1.32); A2.legend(frameon=False, fontsize=5.0, loc='upper right')
A1.text(-0.34, 1.10, 'a', transform=A1.transAxes, fontsize=8, va='bottom')
despine(A1); despine(A2)
tau = pd.read_csv(f'{R}/interface_gtex_tau.csv', index_col=0)['tau'].dropna()
xj = np.arange(9); rj = np.random.default_rng(0)
for i, o in enumerate(ORGS):
    tt = [g for g in M.index if M.loc[g, o] == 1 and g in tau.index]
    v = tau.loc[tt].values
    B.scatter(np.full(len(v), i) + rj.uniform(-0.18, 0.18, len(v)), v, s=6, color=CB[2], alpha=0.5, zorder=2, linewidths=0)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    B.plot([i, i], [q1, q3], color='#9A9A9A', lw=3.2, alpha=0.45, solid_capstyle='round', zorder=3)
    B.scatter([i], [med], s=26, color=CB[1], zorder=5, linewidths=0)
B.set_xticks(xj); B.set_xticklabels(labs, rotation=45, ha='right', fontsize=5.6)
B.set_ylabel('GTEx tau of interface genes'); B.set_ylim(0.70, 1.03)
B.text(0.98, 0.04, 'grey bar = IQR, dot = median', transform=B.transAxes, ha='right', fontsize=5.0, color='#888888')
core = [g for g in M[M['n_organs'] >= 4].index if g in set(AGE['gene'])]
cv = AGE.set_index('gene').loc[core, ['n_sig_pos', 'mean_rho']].sort_values('n_sig_pos')
bg_med = AGE[~AGE['is_iface']]['n_sig_pos'].median()
if_med = AGE[AGE['is_iface']]['n_sig_pos'].median()
yy = np.arange(len(cv))
C.axvspan(0, bg_med, color='#EFEFEF', zorder=0)
C.axvspan(bg_med, if_med, color='#E7F1F8', zorder=0)
C.axvline(bg_med, color='#AAAAAA', lw=0.8, ls=':', zorder=1)
C.axvline(if_med, color='#4A90C2', lw=0.9, ls='--', zorder=1)
C.hlines(yy, 0, cv['n_sig_pos'].values, color='#CCCCCC', lw=1.4, zorder=2)
cc = ['#D55E00' if r > 0 else '#3B7EA1' for r in cv['mean_rho']]
C.scatter(cv['n_sig_pos'], yy, s=52, c=cc, zorder=4, edgecolors='white', linewidths=0.7)
C.text(bg_med / 2, -1.0, 'background\nmedian 3', ha='center', va='top', fontsize=4.8, color='#999999')
C.text((bg_med + if_med) / 2, -1.0, 'interface\nmedian 10', ha='center', va='top', fontsize=4.8, color='#4A90C2')
C.set_yticks(yy); C.set_yticklabels(cv.index, fontsize=5.6)
C.set_ylim(-2.4, len(cv) - 1 + 2.6); C.set_xlim(0, max(cv['n_sig_pos'].max(), if_med) * 1.14)
C.set_xlabel('# tissues with age↑ (of 49)')
C.legend(handles=[Line2D([], [], marker='o', ls='', mfc='#D55E00', mec='none', ms=5, label='mean rho > 0'),
                  Line2D([], [], marker='o', ls='', mfc='#3B7EA1', mec='none', ms=5, label='mean rho < 0')],
         frameon=False, fontsize=4.8, loc='upper right')

ven = pd.read_csv(f'{R}/valid_enrichment.csv'); ven = ven[ven['method'] == 'threshold'].set_index('organ').reindex(ORGS)
_b = pd.read_csv(f'{R}/valid_enrichment.csv').query("organ=='brain28146' and method=='top5%'")['fold']
if len(_b): ven.loc['brain', 'fold'] = float(_b.iloc[0])
D.plot([0, 6.0], [0, 6.0], color='#BBBBBB', lw=0.8, ls='--', zorder=1)
D.scatter(en['fold'], ven['fold'], s=28, color=CB[1], zorder=3, linewidths=0)
bi4 = ORGS.index('brain')
D.scatter([en['fold'].iloc[bi4]], [ven['fold'].iloc[bi4]], s=58, color=CB[2], marker='D', zorder=5, linewidths=0)
for i, o in enumerate(ORGS):
    D.annotate(LBL[o], (en['fold'].iloc[i], ven['fold'].iloc[i]), fontsize=5.2, xytext=(3, 2), textcoords='offset points')
D.set_xlabel('discovery cohort fold'); D.set_ylabel('independent cohort fold')
D.set_xlim(0, 6.0); D.set_ylim(0, 6.0)
for ax, L in zip([B, C, D], 'bcd'): despine(ax); plab(ax, L)
save(fig, 'figure4_validation')
print('AD independent fold =', float(ven.loc['brain', 'fold']))
