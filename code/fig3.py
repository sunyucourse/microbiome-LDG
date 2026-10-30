#!/usr/bin/env python3
"""Figure 3 (manuscript): H2 endotherm vs ectotherm host LDG.
Layout: row 1 = a) gut curves, b) skin curves; row 2 = c) forest plot (centered).
Inputs: gam_h2_*_{curves,strength}.csv
Outputs: /workspace/figures_ms/fig3_host_thermo.{svg,png}
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams.update({'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13,
                            'xtick.labelsize': 11.5, 'ytick.labelsize': 11.5,
                            'legend.fontsize': 11.5})
import matplotlib.pyplot as plt

PROC = '/workspace/data/processed'
C_ENDO = '#D55E00'
C_ECTO = '#0072B2'

def load(prefix):
    c = pd.read_csv(f'{PROC}/{prefix}_curves.csv')
    s = pd.read_csv(f'{PROC}/{prefix}_strength.csv')
    return c, s

fig = plt.figure(figsize=(11, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.35, wspace=0.3)

for pi, (site, title) in enumerate([('gut', 'a  Gut'), ('skin', 'b  Skin / surface')]):
    ax = fig.add_subplot(gs[0, pi])
    curves, strength = load(f'gam_h2_{site}')
    for thermo, col, lab in [('endotherm', C_ENDO, 'Endotherm hosts'),
                             ('ectotherm', C_ECTO, 'Ectotherm hosts')]:
        c = curves[curves['habitat'] == f'{site}_{thermo}']
        if len(c) == 0:
            continue
        ax.fill_between(c['lat_abs'], c['fit'] - 1.96 * c['se'], c['fit'] + 1.96 * c['se'],
                        color=col, alpha=0.2, lw=0)
        ax.plot(c['lat_abs'], c['fit'], color=col, lw=2.2, label=lab)
    ax.set_xlabel('|Latitude| (degrees)')
    if pi == 0:
        ax.set_ylabel('log richness (GAM fit)')
    ax.set_title(title, loc='left', fontweight='bold')
    ax.legend(frameon=False, loc='best')

# forest plot, centered in row 2
ax = fig.add_subplot(gs[1, :])
rows = []
for site in ['gut', 'skin']:
    for suffix, lab in [('', 'all hosts'), ('_wildonly', 'wild only')]:
        fn = f'{PROC}/gam_h2_{site}{suffix}_strength.csv'
        if not os.path.exists(fn):
            continue
        s = pd.read_csv(fn)
        if 'n_grids' not in s.columns:
            continue
        for _, r in s.iterrows():
            thermo = r['habitat'].split('_')[-1]
            rows.append(dict(site=site, thermo=thermo, lab=lab,
                             ratio=r['eq_vs_highlat_ratio'],
                             lo=r['ratio_lo'], hi=r['ratio_hi'],
                             n=r['n_grids']))
f = pd.DataFrame(rows)
f = f.iloc[::-1].reset_index(drop=True)
for i, r in f.iterrows():
    col = C_ENDO if r['thermo'] == 'endotherm' else C_ECTO
    if np.isfinite(r['ratio']):
        ax.plot([r['lo'], r['hi']], [i, i], color=col, lw=1.8)
        ax.plot(r['ratio'], i, 'o', color=col, ms=8)
    ax.text(0.16, i + 0.34, f"{r['site']} {r['thermo']} ({r['lab']}, n={int(r['n'])})",
            ha='left', va='bottom', fontsize=11)
ax.axvline(1, color='grey', lw=0.8, ls='--')
ax.set_xscale('log')
ax.set_xlim(0.15, 15)
ax.set_xticks([0.2, 0.5, 1, 2, 5, 10])
ax.set_xticklabels(['0.2', '0.5', '1', '2', '5', '10'])
ax.set_ylim(-0.6, len(f) - 0.1)
ax.set_yticks([])
ax.set_xlabel('Equator (5°) / high-latitude (55°) richness ratio')
ax.set_title('c  Gradient strength by host thermoregulation', loc='left', fontweight='bold')
ax.spines[['left']].set_visible(False)

fig.savefig('/workspace/figures_ms/fig3_host_thermo.svg', bbox_inches='tight')
fig.savefig('/workspace/figures_ms/fig3_host_thermo.png', bbox_inches='tight', dpi=200)
print('fig3 saved')
