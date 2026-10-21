"""Build the lineage-split IDO1 trajectory figures (x2, x4) as vector SVGs
with editable text, reproducing the canonical style from committed data.
No suptitles (captions carry the stats).

Inputs:
  hub_spot_candidate_genes_with_pt.csv  (37,230 Hub/TLS spots; IDO1, pseudotime_l1/l2, tls_stage, cancer_type)
  ido1_pseudotime_values.csv            (116,942 scRNA cells; IDO1, pseudotime)

Outputs (to metabolism_target/ido1_report_figures/):
  fig_x2_ridgeline.svg / .png   L1 ridgeline + cancer x stage heatmap
  fig_x4_loess.svg / .png       dual trajectory: scRNA fraction + spatial L1/L2 fractions
Prints all computed stats.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import stats

BASE = 'data/hub_tls_paper_outline/metabolism_target'
OUT = f'{BASE}/ido1_report_figures'
plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.bbox': 'tight'})

spots = pd.read_csv(f'{BASE}/hub_spot_candidate_genes_with_pt.csv')
ptc = pd.read_csv(f'{BASE}/ido1_pseudotime_values.csv')

BINS = [0, 7, 11, 13, 17, 20, 24, 32, 38]
BIN_LABS = [f'{BINS[i]}-{BINS[i+1]}' for i in range(len(BINS)-1)]

def frac_ci(k, n, z=1.96):
    p = k / n
    denom = 1 + z**2 / n
    c = (p + z**2/(2*n)) / denom
    h = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
    return (c-h)*100, (c+h)*100

# ===================================================================
# Fig x4: dual trajectory (canonical layout: left scRNA, right spatial L1/L2)
# ===================================================================
C_SCRNA = '#b83a4b'; C_L1 = '#2b6cb0'; C_L2 = '#c07a2b'

def frac_curve(ax, x, y, color, label, linestyle='-', nbin=25):
    df = pd.DataFrame({'x': x, 'y': (y > 0).astype(float)})
    df['bin'] = pd.qcut(df['x'], nbin, labels=False, duplicates='drop')
    g = df.groupby('bin').agg(xc=('x','mean'), n=('y','size'), k=('y','sum'))
    lo, hi = frac_ci(g['k'].values, g['n'].values)
    ax.fill_between(g['xc'], lo, hi, color=color, alpha=.15, lw=0)
    ax.plot(g['xc'], g['k']/g['n']*100, linestyle, color=color, lw=2.2,
            marker='o', ms=3.5, label=label)

l1 = spots.dropna(subset=['pseudotime_l1'])
l2 = spots.dropna(subset=['pseudotime_l2'])
rho_s = stats.spearmanr(ptc['pseudotime'], ptc['IDO1'])[0]
rho1 = stats.spearmanr(l1['IDO1'], l1['pseudotime_l1'])[0]
rho2 = stats.spearmanr(l2['IDO1'], l2['pseudotime_l2'])[0]

fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
frac_curve(axes[0], ptc['pseudotime'].values, ptc['IDO1'].values, C_SCRNA,
           f'scRNA cells (n = {len(ptc):,})')
axes[0].set_xlabel('scRNA-seq pseudotime (cell trajectory)')
axes[0].set_ylabel('IDO1+ cells/spots (%)')
axes[0].annotate(f'ρ = +{rho_s:.3f}', xy=(.03, .97), xycoords='axes fraction',
                 fontsize=8.5, va='top', color=C_SCRNA)
frac_curve(axes[1], l1['pseudotime_l1'].values, l1['IDO1'].values, C_L1,
           f'Lineage 1 → SFL-TLS endpoint (ρ = {rho1:.3f})', '-')
frac_curve(axes[1], l2['pseudotime_l2'].values, l2['IDO1'].values, C_L2,
           f'Lineage 2 → E-TLS endpoint (ρ = {rho2:.3f})', '--')
axes[1].set_xlabel('Spatial pseudotime (per lineage)')
axes[1].legend(frameon=False, fontsize=7.5, loc='upper right')
plt.tight_layout()
plt.savefig(f'{OUT}/fig_x4_loess.svg'); plt.savefig(f'{OUT}/fig_x4_loess.png', dpi=300)
plt.close()
print(f'[x4] scRNA rho=+{rho_s:.4f} (n={len(ptc)}); L1 rho={rho1:.4f} (n={len(l1)}); L2 rho={rho2:.4f} (n={len(l2)})')

# ===================================================================
# Fig x2: L1 ridgeline (custom bins, IDO1+% annotated) + cancer x stage heatmap
# ===================================================================
l1b = l1.copy()
l1b['bin'] = pd.cut(l1b['pseudotime_l1'], bins=BINS, labels=BIN_LABS, include_lowest=True)
l1b = l1b.dropna(subset=['bin'])
fracs = l1b.groupby('bin', observed=True)['IDO1'].apply(lambda s: (s > 0).mean() * 100)

fig, (ax_r, ax_h) = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={'width_ratios': [1, 1.35]})
# ridgeline (vertical kde per bin) with a GLOBAL x scale
from scipy.stats import gaussian_kde
nb = len(BIN_LABS)
xmax = float(np.percentile(l1b['IDO1'], 99))
xs = np.linspace(0, max(xmax, 0.5), 200)
for i, lab in enumerate(BIN_LABS):
    v = l1b.loc[l1b['bin'] == lab, 'IDO1'].values
    y0 = nb - 1 - i
    if len(v) > 10 and v.std() > 0:
        kde = gaussian_kde(v)(xs)
        kde = kde / kde.max() * 0.85
        ax_r.fill_betweenx(y0 + kde, 0, xs, color=cm.plasma(0.15 + 0.7 * i / (nb - 1)), alpha=.9, lw=0)
        ax_r.plot(xs, y0 + kde, color=cm.plasma(0.15 + 0.7 * i / (nb - 1)), lw=.6)
    ax_r.text(xmax * 1.03, y0, f'{lab}\nIDO1+ {fracs[lab]:.1f}%',
              va='center', fontsize=7, color=cm.plasma(0.15 + 0.7 * i / (nb - 1)))
ax_r.set_xlim(0, xmax * 1.42)
ax_r.set_ylim(-.6, nb - .4)
ax_r.set_yticks([])
ax_r.set_xlabel('IDO1 expression')
ax_r.set_ylabel('Spatial pseudotime bin — Lineage 1 (ends in SFL-TLS)')

# heatmap: mean IDO1 + IDO1+% per cancer x stage (rows ordered by E-TLS, canonical order)
stages = ['E-TLS', 'PFL-TLS', 'SFL-TLS']
piv_m = spots.pivot_table(index='cancer_type', columns='tls_stage', values='IDO1', aggfunc='mean')[stages]
order = piv_m['E-TLS'].sort_values(ascending=False).index
piv_m = piv_m.reindex(order)
piv_f = spots.assign(p=spots['IDO1'] > 0).pivot_table(index='cancer_type', columns='tls_stage', values='p', aggfunc='mean').reindex(order)[stages] * 100
im = ax_h.imshow(piv_m.values, cmap='magma', aspect='auto')
ax_h.set_xticks(range(3)); ax_h.set_xticklabels(stages, fontsize=8)
ax_h.set_yticks(range(len(order))); ax_h.set_yticklabels(order, fontsize=8)
for i in range(len(order)):
    for j in range(3):
        v = piv_m.values[i, j]
        if not np.isnan(v):
            ax_h.text(j, i, f'{v:.2f}\n({piv_f.values[i, j]:.0f}%)', ha='center', va='center',
                      fontsize=6.5, color='white' if v < piv_m.values[~np.isnan(piv_m.values)].mean() else 'black')
fig.colorbar(im, ax=ax_h, label='Mean IDO1', shrink=.8)
plt.tight_layout()
plt.savefig(f'{OUT}/fig_x2_ridgeline.svg'); plt.savefig(f'{OUT}/fig_x2_ridgeline.png', dpi=300)
plt.close()

print('[x2] L1 bins (n={}): IDO1+% ='.format(len(l1b)), fracs.round(1).tolist())
print('[x2] L1 rho={:.4f} (n={})'.format(rho1, len(l1)))
