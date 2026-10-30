#!/usr/bin/env python3
"""Figure 2: per-habitat LDG curves ranked by gradient strength + H1 summary.
Inputs: gam_all_curves.csv, gam_all_strength.csv (from fit_ldg_gam.R)
Outputs: /mnt/results/figures/fig2_ldg_curves.{svg,png}
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams.update({'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13,
                            'xtick.labelsize': 11.5, 'ytick.labelsize': 11.5})
import matplotlib.pyplot as plt

sys.path.insert(0, '/workspace/ldg')
from figstyle import HABITAT_COLORS, HABITAT_LABELS, EXPOSURE

PROC = '/workspace/data/processed'
curves = pd.read_csv(f'{PROC}/gam_all_curves.csv')
strength = pd.read_csv(f'{PROC}/gam_all_strength.csv')
strength['exposure'] = strength['habitat'].map(EXPOSURE)

# order habitats by lat2 coefficient (latitudinal-structure strength; negative = peak/decline structure)
strength = strength.sort_values('lat2_coef', ascending=True)
order = strength['habitat'].tolist()

fig = plt.figure(figsize=(11, 9))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.05)

# --- panel a: GAM curves, vertically offset in strength order ---
axa = fig.add_subplot(gs[0, 0])
offset_step = 1.0
for i, hab in enumerate(order):
    c = curves[curves['habitat'] == hab]
    if len(c) == 0:
        continue
    off = -i * offset_step
    y = c['fit'] - c['fit'].median() + off
    axa.fill_between(c['lat_abs'], y - 1.96 * c['se'], y + 1.96 * c['se'],
                     color=HABITAT_COLORS[hab], alpha=0.25, lw=0)
    axa.plot(c['lat_abs'], y, color=HABITAT_COLORS[hab], lw=1.8)
    axa.text(-1.5, off, HABITAT_LABELS[hab], ha='right', va='center', fontsize=11,
             color=HABITAT_COLORS[hab])
axa.set_xlim(0, 70)
axa.set_ylim(-len(order) * offset_step, offset_step)
axa.set_yticks([])
axa.set_xlabel('|Latitude| (degrees)')
axa.set_ylabel('log richness (centered per habitat)', labelpad=145)
axa.yaxis.set_label_coords(-0.32, 0.5)
axa.set_title('a  LDG curves ranked by strength', loc='left', fontweight='bold')
axa.spines[['left']].set_visible(False)

# --- panel b: lat2 coefficient (gradient strength) with SE-based CI ---
axb = fig.add_subplot(gs[0, 1])
yy = -np.arange(len(order)) * offset_step
for i, (_, r) in enumerate(strength.iterrows()):
    hab = r['habitat']
    lo, hi = r['lat2_coef'] - 1.96 * r['lat2_se'], r['lat2_coef'] + 1.96 * r['lat2_se']
    axb.plot([lo, hi], [yy[i], yy[i]], color=HABITAT_COLORS[hab], lw=1.5)
    axb.plot(r['lat2_coef'], yy[i], 'o', color=HABITAT_COLORS[hab], ms=6)
axb.axvline(0, color='grey', lw=0.8, ls='--')
axb.set_ylim(-len(order) * offset_step, offset_step)
axb.set_yticks([])
axb.set_xlabel('Latitudinal-structure coefficient (lat², quadratic GLM)\nnegative = mid-latitude peak / decline from peak')
axb.set_title('b  Gradient strength', loc='left', fontweight='bold')
axb.spines[['left']].set_visible(False)

fig.savefig('/workspace/figures_ms/fig2_ldg_curves.svg', bbox_inches='tight')
fig.savefig('/workspace/figures_ms/fig2_ldg_curves.png', bbox_inches='tight', dpi=200)
print('fig2 saved')
