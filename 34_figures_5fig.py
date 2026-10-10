import json, os, math
from collections import Counter
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Circle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import networkx as nx

R = 'out'; OUT = 'out/figures'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'], 'font.size': 7,
                     'axes.linewidth': 0.6, 'axes.labelsize': 6.2, 'xtick.labelsize': 6.2, 'ytick.labelsize': 6.2,
                     'legend.fontsize': 6.2, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5})
CB = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#999999']
ORGS = ['kidney', 'liver', 'heart', 'lung', 'colon', 'skin', 'artery', 'esophagus', 'brain']
LBL = {'kidney': 'DKD', 'liver': 'Cirrhosis', 'heart': 'DCM', 'lung': 'IPF', 'colon': 'IBD',
       'skin': 'Psoriasis', 'artery': 'Athero', 'esophagus': 'EAC', 'brain': 'AD'}
labs = [LBL[o] for o in ORGS]; x = np.arange(9)
SIGORD = ['SenCat', 'SenMayo', 'Fridman', 'ReactomeSASP', 'Reactome', 'GO0090398']

en = pd.read_csv(f'{R}/interface_enrichment.csv').set_index('organ').loc[ORGS]
M = pd.read_csv(f'{R}/interface_matrix.csv', index_col=0)
sc = pd.read_csv(f'{R}/interface_score_per_organ.csv').set_index('organ').loc[ORGS]
null = json.load(open(f'{R}/crossorgan_null.json'))
ms = pd.read_csv(f'{R}/multisig_enrichment.csv')
AGE = pd.read_csv(f'{R}/gtex_age_summary.csv')
path_s = pd.read_csv(f'{R}/pathway_shared_ge3organs.csv'); path_o = pd.read_csv(f'{R}/pathway_organ_specific.csv')
sub = pd.read_csv(f'{R}/subtype_consistency_summary.csv').set_index('organ')
RD = pd.read_csv(f'{R}/brain_region_sencat_ci.csv')

def despine(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
def plab(ax, L):
    ax.text(-0.16, 1.05, L, transform=ax.transAxes, fontsize=8, va='bottom', ha='left')
def setx(ax, size):
    for t in ax.get_xticklabels(): t.set_fontsize(size); t.set_fontweight('normal')
def save(fig, name):
    fig.savefig(f'{OUT}/{name}.pdf', dpi=600)
    fig.savefig(f'{OUT}/{name}.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig); print('saved', name)

fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.34, left=0.09, right=0.975, top=0.94, bottom=0.13)
A, B, C = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0)]]
A.bar(x, en['fold'], color=[CB[i % 5] for i in range(9)], width=0.7)
A.axhline(1, color='gray', lw=0.6, ls='--')
for i, (f, o) in enumerate(zip(en['fold'], en['overlap'])): A.text(i, f + 0.12, '%d' % o, ha='center', fontsize=5.2)
A.set_xticks(x); A.set_xticklabels(labs, rotation=45, ha='right', fontsize=5.8)
A.set_ylabel('enrichment (obs/exp)'); A.set_ylim(0, en['fold'].max() * 1.3)
ks = [2, 3, 4, 5]; obs = [null['obs'][str(k)] for k in ks]; nul = [null['nullmean'][str(k)] for k in ks]
B.bar(ks, obs, width=0.55, color=CB[0], label='observed'); B.scatter(ks, nul, marker='D', s=20, color=CB[6], zorder=5, label='null (random)')
for k, o in zip(ks, obs): B.text(k, o * 1.25, str(o), ha='center', fontsize=6)
B.set_yscale('symlog', linthresh=1, linscale=0.7, subs=[2, 5]); B.set_ylim(0, 400)
B.set_xticks(ks); B.set_xticklabels(['≥2', '≥3', '≥4', '≥5']); B.set_xlabel('shared by ≥k organs'); B.set_ylabel('# interface genes')
B.text(0.98, 0.62, 'all p<0.001', transform=B.transAxes, ha='right', fontsize=6.2)
B.legend(frameon=False, loc='upper right')
C.bar(x, sc['case_med'], width=0.62, color=CB[1], label='disease'); C.bar(x, sc['ctrl_med'], width=0.62, color=CB[0], alpha=0.55, label='control')
C.axhline(0, color='black', lw=0.5); C.set_xticks(x); C.set_xticklabels(labs, rotation=45, ha='right', fontsize=5.8)
C.set_ylabel('interface module score'); C.legend(frameon=False, loc='upper right')
gsd = gs[1, 1].subgridspec(2, 1, height_ratios=[1.5, 2.0], hspace=0.10)
axT = fig.add_subplot(gsd[0]); axBm = fig.add_subplot(gsd[1], sharex=axT)
combos = Counter(frozenset(o for o in ORGS if M.loc[g, o] == 1) for g in M.index)
multi = sorted([(c, n) for c, n in combos.items() if len(c) >= 2], key=lambda t: -t[1])[:16]
sing = sum(n for c, n in combos.items() if len(c) == 1)
xs = np.arange(len(multi))
axT.bar(xs, [n for _, n in multi], color=CB[0], width=0.62)
for i, (_, n) in enumerate(multi): axT.text(i, n + 0.15, str(n), ha='center', fontsize=6.2)
axT.set_ylabel('# interface genes'); axT.set_ylim(0, max(n for _, n in multi) * 1.25)
axT.text(0.98, 0.80, '%d genes are\norgan-specific' % sing, transform=axT.transAxes, ha='right', fontsize=6, color='#666666')
despine(axT); axT.tick_params(bottom=False, labelbottom=False)
for j, o in enumerate(ORGS): axBm.scatter(xs, np.full(len(xs), j), s=13, color='#E5E5E5', zorder=1)
for i, (c, n) in enumerate(multi):
    ys = [ORGS.index(o) for o in c]
    for y in ys: axBm.scatter(i, y, s=20, color=CB[1], zorder=3)
    if len(ys) > 1: axBm.plot([i, i], [min(ys), max(ys)], color=CB[1], lw=0.9, zorder=2)
axBm.set_yticks(range(len(ORGS))); axBm.set_yticklabels(labs, fontsize=5)
axBm.set_ylim(-0.6, len(ORGS) - 0.4)
axBm.set_xticks(xs); axBm.set_xticklabels([])
axBm.set_xlabel('organ combinations sharing a gene', fontsize=6.2)
for s in ['top', 'right', 'left']: axBm.spines[s].set_visible(False)
axBm.tick_params(axis='y', length=0)
for ax, L in zip([A, B, C], 'abc'): despine(ax); plab(ax, L)
plab(axT, 'd')
setx(A, 6.2); setx(B, 6.2); setx(C, 6.2); setx(axBm, 5)
save(fig, 'figure1_discovery')

subM = M[M['n_organs'] >= 3]
G = nx.Graph()
for o in ORGS: G.add_node(o, kind='organ')
for g in subM.index: G.add_node(g, kind='gene', n=int(subM.loc[g, 'n_organs']))
for g in subM.index:
    for o in ORGS:
        if subM.loc[g, o] == 1: G.add_edge(o, g)
ang = {o: 2 * math.pi * j / len(ORGS) + math.pi / 2 for j, o in enumerate(ORGS)}
pos = {o: (math.cos(ang[o]), math.sin(ang[o])) for o in ORGS}
seen = {}
for g in [n for n in G if G.nodes[n]['kind'] == 'gene']:
    Aa = [ang[o] for o in ORGS if subM.loc[g, o] == 1]
    a = math.atan2(np.mean([math.sin(t) for t in Aa]), np.mean([math.cos(t) for t in Aa]))
    n = G.nodes[g]['n']; r = 0.90 - 0.095 * n
    k = seen.get(round(a, 2), 0); seen[round(a, 2)] = k + 1
    r += 0.05 * (k % 4); a += 0.07 * (k // 4)
    pos[g] = (r * math.cos(a), r * math.sin(a))
fig, ax = plt.subplots(figsize=(6.0, 6.4))
fig.subplots_adjust(left=0.04, right=0.96, top=0.99, bottom=0.115)
nx.draw_networkx_edges(G, pos, ax=ax, edge_color='#D9D9D9', width=0.4, alpha=0.55)
genes = [n for n in G if G.nodes[n]['kind'] == 'gene']
nx.draw_networkx_nodes(G, pos, nodelist=genes, ax=ax, node_size=[26 * G.nodes[n]['n'] for n in genes],
                       node_color=[CB[1] if G.nodes[n]['n'] >= 5 else (CB[0] if G.nodes[n]['n'] == 4 else '#9ECAE1') for n in genes],
                       alpha=0.92, linewidths=0)
R_ORG = 0.17
for o in ORGS:
    ax.add_patch(Circle(pos[o], R_ORG, facecolor='white', edgecolor=CB[2], lw=1.5, linestyle='-', zorder=6))
    ax.text(pos[o][0], pos[o][1], LBL[o], ha='center', va='center', fontsize=8.5, fontweight='bold', zorder=7)
for g in genes:
    if G.nodes[g]['n'] >= 5:
        ax.text(pos[g][0], pos[g][1], g, ha='center', va='center', fontsize=7.0, fontweight='bold', color='black',
                path_effects=[pe.withStroke(linewidth=2.4, foreground='white')])
ax.set_xlim(-1.22, 1.22); ax.set_ylim(-1.22, 1.22); ax.set_aspect('equal'); ax.axis('off')
ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=CB[1], mec='none', ms=8, label='shared by ≥5 organs (core)'),
                   Line2D([], [], marker='o', ls='', mfc=CB[0], mec='none', ms=7, label='shared by 4 organs'),
                   Line2D([], [], marker='o', ls='', mfc='#9ECAE1', mec='none', ms=6, label='shared by 3 organs')],
          frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.005), ncol=3, fontsize=8,
          handletextpad=0.4, columnspacing=1.4)
save(fig, 'figure2_network')

fig = plt.figure(figsize=(7.2, 4.5))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.36, left=0.15, right=0.975, top=0.94, bottom=0.115)
A = fig.add_subplot(gs[0, 0]); B = fig.add_subplot(gs[0, 1]); C = fig.add_subplot(gs[1, 0]); D = fig.add_subplot(gs[1, 1])
tops = pd.concat([path_s.sort_values('p').head(5).assign(cat='shared'), path_o.sort_values('p').head(5).assign(cat='specific')])
tops['nlp'] = -np.log10(tops['p']); tops = tops.sort_values('nlp')
A.barh(range(len(tops)), tops['nlp'], color=[CB[1] if c == 'shared' else CB[0] for c in tops['cat']], height=0.72)
A.set_yticks(range(len(tops)))
A.set_yticklabels([(t if len(t) <= 24 else t[:24] + '…' if t[24] == ' ' else t[:24].rsplit(' ', 1)[0] + '…') for t in tops['term']], fontsize=5.4)
A.set_xlabel('-log10 p')
A.legend(handles=[Patch(facecolor=CB[1], label='shared by ≥3 organs'), Patch(facecolor=CB[0], label='organ-specific')],
         frameon=False, loc='lower right', fontsize=5.6)
cnt = ms[ms['p'] < 0.05].groupby('sig')['organ'].nunique().reindex(SIGORD).fillna(0)
B.bar(range(len(SIGORD)), cnt.values, color=[CB[i % len(CB)] for i in range(len(SIGORD))], width=0.6)
for i, v in enumerate(cnt.values): B.text(i, v + 0.1, str(int(v)), ha='center', fontsize=6)
B.set_xticks(range(len(SIGORD))); B.set_xticklabels(SIGORD, rotation=45, ha='right', fontsize=6.2)
B.set_ylabel('# organs enriched (p<0.05)'); B.set_ylim(0, 9)
C.scatter(en['n_sig'], en['fold'], s=30, color=CB[2], zorder=5)
for _, r in en.iterrows():
    lab = LBL[r.name]
    if lab == 'Athero':
        C.annotate(lab, (r['n_sig'], r['fold']), fontsize=5.4, xytext=(4, 0), textcoords='offset points', ha='left', va='center')
    else:
        C.annotate(lab, (r['n_sig'], r['fold']), fontsize=5.4, xytext=(-4, 3), textcoords='offset points', ha='right')
C.set_xscale('log'); C.set_xlabel('# disease-signature genes (log)'); C.set_ylabel('enrichment fold')
C.axhline(1, color='gray', lw=0.5, ls=':')
tau = pd.read_csv(f'{R}/interface_gtex_tau.csv', index_col=0)['tau']
tm = M[['n_organs']].join(tau, how='inner').dropna()
tm['grp'] = tm['n_organs'].clip(upper=3)
data = [tm.loc[tm['grp'] == g, 'tau'].values for g in (1, 2, 3)]
vp = D.violinplot(data, positions=range(3), widths=0.78, showextrema=False, showmedians=True)
for pc, c in zip(vp['bodies'], [CB[0], CB[2], CB[3]]):
    pc.set_facecolor(c); pc.set_alpha(0.65); pc.set_edgecolor('#555555'); pc.set_linewidth(0.6)
vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.0)
from scipy.stats import spearmanr
rho, pv = spearmanr(tm['n_organs'], tm['tau'])
D.text(0.97, 0.96, 'Spearman ρ=%.2f, p=%.1g' % (rho, pv), transform=D.transAxes, ha='right', va='top', fontsize=5.6, color='#666666')
D.set_xticks(range(3)); D.set_xticklabels(['%s\n(n=%d)' % (g, len(d)) for g, d in zip(['1', '2', '≥3'], data)])
D.set_xlabel('# organs sharing the gene')
D.set_ylabel('GTEx tau (higher = tissue-restricted)'); D.set_ylim(0.70, 1.06)
for ax, L in zip([A, B, C, D], 'abcd'): despine(ax); plab(ax, L)
save(fig, 'figure3_mechanism')

fig = plt.figure(figsize=(7.2, 4.4))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.32, left=0.09, right=0.975, top=0.90, bottom=0.115)
gsa = gs[0, 0].subgridspec(1, 2, wspace=0.42)
A1 = fig.add_subplot(gsa[0]); A2 = fig.add_subplot(gsa[1])
B = fig.add_subplot(gs[0, 1]); C = fig.add_subplot(gs[1, 0]); D = fig.add_subplot(gs[1, 1])
from scipy.stats import gaussian_kde
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
C.set_ylim(-2.4, 6.4)
C.set_xlim(0, max(cv['n_sig_pos'].max(), if_med) * 1.14)
C.set_xlabel('# tissues with age↑ (of 49)')
C.legend(handles=[Line2D([], [], marker='o', ls='', mfc='#D55E00', mec='none', ms=5, label='mean rho > 0'),
                  Line2D([], [], marker='o', ls='', mfc='#3B7EA1', mec='none', ms=5, label='mean rho < 0')],
         frameon=False, fontsize=4.8, loc='upper right')
ven = pd.read_csv(f'{R}/valid_enrichment.csv'); ven = ven[ven['method'] == 'threshold'].set_index('organ').reindex(ORGS)
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

fig = plt.figure(figsize=(7.2, 3.0))
gs = fig.add_gridspec(1, 2, wspace=0.34, left=0.085, right=0.985, top=0.85, bottom=0.22)
axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1])
FO = ms.pivot_table(index='organ', columns='sig', values='fold').reindex(ORGS)[SIGORD]
ad = FO.loc['brain']; oth = FO.drop(index='brain').median()
xs = np.arange(len(SIGORD))
axa.bar(xs - 0.2, oth.values, width=0.4, color='#BBBBBB', label='other 8 organs (median)')
axa.bar(xs + 0.2, ad.values, width=0.4, color=CB[1], label='brain (AD)')
axa.axhline(1, color='gray', lw=0.6, ls='--')
axa.set_xticks(xs); axa.set_xticklabels(SIGORD, rotation=45, ha='right', fontsize=6.2)
axa.set_ylabel('enrichment fold'); axa.legend(frameon=False, fontsize=4.8, loc='upper right')
axa.set_ylim(0, max(ad.max(), oth.max()) * 1.4)
axb.errorbar(RD['fold'], range(len(RD)), xerr=[RD['fold'] - RD['lo'], RD['hi'] - RD['fold']], fmt='o', color=CB[0], ecolor=CB[0], elinewidth=0.9, capsize=2, ms=3.4)
for i, (hi, cap) in enumerate(zip(RD['hi'], RD['capped'])):
    if cap: axb.annotate('', xy=(hi * 1.06, i), xytext=(hi, i), arrowprops=dict(arrowstyle='->', color=CB[0], lw=0.9))
axb.axvline(1, color='gray', lw=0.7, ls='--')
axb.set_yticks(range(len(RD))); axb.set_yticklabels([r.replace('_', ' ') for r in RD['region']], fontsize=6.2)
axb.set_xlabel('SenCat fold (per brain region)'); axb.set_xlim(-0.12, 4.15)
for i, p in enumerate(RD['p']): axb.text(3.5, i, 'p=%.2f' % p, va='center', ha='left', fontsize=4.8, color='#666666')
for ax, L in zip([axa, axb], 'ab'): despine(ax); plab(ax, L)
save(fig, 'figure5_brain')
print('DONE ->', OUT)
