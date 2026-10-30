#!/usr/bin/env python3
"""Figure 4: H3 engineered habitats decoupled from latitude.
Panels: a) partial R2 of |lat| vs human-population/GDP/sanitation predictors;
        b) same-latitude contrast: engineered vs natural richness log-ratio
           across latitude.
Inputs: h3_variance_decomposition.csv, h3_contrast.csv
Outputs: /mnt/results/figures/fig4_engineered_decoupling.{svg,png}
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams.update({'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13, 'xtick.labelsize': 11.5, 'ytick.labelsize': 11.5, 'legend.fontsize': 11.5})
import matplotlib.pyplot as plt

sys.path.insert(0, '/workspace/ldg')
from figstyle import HABITAT_COLORS, HABITAT_LABELS

PROC = '/workspace/data/processed'
vd = pd.read_csv(f'{PROC}/h3_variance_decomposition.csv')
con = pd.read_csv(f'{PROC}/h3_contrast.csv')

fig = plt.figure(figsize=(11, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.3)

# --- panel a: partial R2 bars ---
axa = fig.add_subplot(gs[0, 0])
preds = [('lat_abs', '|Latitude|', '#4C72B0'),
         ('log_pop', 'log population', '#DD8452'),
         ('log_gdp', 'log GDP/capita', '#55A868'),
         ('sanit', 'basic sanitation', '#C44E52')]
habs = vd['habitat'].tolist()
y = np.arange(len(habs))
left = np.zeros(len(habs))
for p, lab, col in preds:
    c = f'partial_R2_{p}'
    vals = vd[c].fillna(0).to_numpy() if c in vd.columns else np.zeros(len(habs))
    axa.barh(y, vals, left=left, color=col, label=lab, height=0.65)
    left += vals
axa.set_yticks(y)
axa.set_yticklabels([HABITAT_LABELS.get(h, h) for h in habs])
axa.invert_yaxis()
axa.set_xlabel('Partial R² (variance in log richness explained)')
axa.legend(frameon=False, fontsize=11, loc='lower right')
axa.set_title('a  What predicts engineered-habitat diversity?', loc='left', fontweight='bold')

# --- panel b: same-latitude contrast ---
axb = fig.add_subplot(gs[0, 1])
for hab in con['habitat'].unique():
    d = con[con['habitat'] == hab]
    axb.scatter(d['lat'], d['log_ratio'], s=14, alpha=0.5, color=HABITAT_COLORS.get(hab, 'grey'),
                label=HABITAT_LABELS.get(hab, hab), lw=0)
    # binned median trend
    bins = np.arange(-60, 75, 10)
    d = d.copy()
    d['bin'] = pd.cut(d['lat'], bins)
    bm = d.groupby('bin', observed=True)['log_ratio'].median()
    xc = [b.mid for b in bm.index]
    axb.plot(xc, bm.values, color=HABITAT_COLORS.get(hab, 'grey'), lw=2)
axb.axhline(0, color='grey', lw=0.8, ls='--')
axb.set_xlabel('Latitude (degrees)')
axb.set_ylabel('log(engineered / natural richness)\nat same latitude')
axb.legend(frameon=False, fontsize=11)
axb.set_title('b  Engineered vs natural diversity, same latitude', loc='left', fontweight='bold')

fig.savefig('/workspace/figures_ms/fig4_engineered_decoupling.svg', bbox_inches='tight')
fig.savefig('/workspace/figures_ms/fig4_engineered_decoupling.png', bbox_inches='tight', dpi=200)
print('fig4 saved')
