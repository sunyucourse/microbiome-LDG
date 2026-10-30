#!/usr/bin/env python3
"""Manuscript figures v3: fig1 with visible land base map; fig2/fig3 in grouped habitat order;
4-type direction typology (classic decline / reverse increase / mid-latitude peak / flat)."""
import sys
sys.path.insert(0, '/workspace/ldg')
from figstyle import *
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

FIG = '/workspace/figures_ms'
gt = pd.read_csv('/mnt/results/tables/grid_table.csv')
p01 = pd.read_csv('/mnt/results/tables/p01_shape_classification.csv')
p02 = pd.read_csv('/mnt/results/tables/p02_gradient_metrics.csv')
p03 = pd.read_csv('/mnt/results/tables/p03_study_re.csv')
p13 = pd.read_csv('/mnt/results/tables/p13_unmapped_sensitivity.csv')
curves = pd.read_csv('/workspace/data/processed/gam_all_curves.csv')

# grouped order: soils / waters / engineered / plants / hosts
ORDER = ['soil_natural','soil_agricultural','freshwater','marine_planktonic','marine_sediment',
         'wastewater','urban_built','plant_rhizosphere','plant_phyllosphere',
         'human_gut','animal_gut','animal_skin']
SHAPE_V3 = {'freshwater':'classic decline','marine_planktonic':'classic decline','wastewater':'classic decline',
            'soil_natural':'reverse increase','human_gut':'reverse increase',
            'plant_rhizosphere':'mid-latitude peak','marine_sediment':'mid-latitude peak (n.s.)',
            'animal_gut':'flat','animal_skin':'flat (n.s.)','plant_phyllosphere':'flat (n.s.)',
            'soil_agricultural':'flat (n.s.)','urban_built':'flat (n.s.)'}

# ---------------- Fig 1: map with land base + sample sizes ----------------
land = gpd.read_file('/workspace/data/external/naturalearth/ne_110m_land.shp')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={'width_ratios':[1.7,1]})
land.plot(ax=ax1, color='#E3E1DC', edgecolor='#B8B5AE', linewidth=0.15, zorder=1)
for hab in ORDER:
    x = gt[gt.habitat == hab]
    ax1.scatter(x.glon + 0.5, x.glat + 0.5, s=5, color=HABITAT_COLORS[hab],
                alpha=0.6, linewidths=0, rasterized=True, zorder=2)
ax1.set_xlim(-180, 180); ax1.set_ylim(-62, 85)
ax1.set_xlabel('Longitude'); ax1.set_ylabel('Latitude')
ax1.set_title('a  Grid cells with >=3 samples (n = 5,507)', loc='left', fontsize=10)

samp = gt.groupby('habitat').n_samples.sum().reindex(ORDER)
y = np.arange(len(ORDER))
ax2.barh(y, samp.values, color=[HABITAT_COLORS[h] for h in ORDER], height=0.7)
ax2.set_yticks(y); ax2.set_yticklabels([HABITAT_LABELS[h] for h in ORDER], fontsize=8.5)
ax2.invert_yaxis()
ax2.set_xscale('log'); ax2.set_xlim(500, 3e5)
ax2.set_xlabel('Samples (log scale)')
ax2.set_title('b  Samples per habitat', loc='left', fontsize=10)
for i, v in enumerate(samp.values):
    ax2.text(v*1.15, i, f'{v:,}', va='center', fontsize=7.5)
fig.tight_layout()
fig.savefig(f'{FIG}/fig1_sample_map.png'); fig.savefig(f'{FIG}/fig1_sample_map.svg')
plt.close(fig)
print('fig1 done')

# ---------------- Fig 2: GAM curves 3x4, grouped order, 4-type labels ----------------
fig, axes = plt.subplots(3, 4, figsize=(11, 7.5), sharex=True)
p01i = p01.set_index('habitat')
for ax, hab in zip(axes.flat, ORDER):
    c = curves[curves.habitat == hab]
    ax.fill_between(c.lat_abs, np.exp(c.fit - 1.96*c.se), np.exp(c.fit + 1.96*c.se),
                    color=HABITAT_COLORS[hab], alpha=0.25, linewidth=0)
    ax.plot(c.lat_abs, np.exp(c.fit), color=HABITAT_COLORS[hab], lw=1.8)
    r = p01i.loc[hab]
    ax.set_yscale('log')
    ax.set_title(HABITAT_LABELS[hab], fontsize=9, loc='left')
    ax.text(0.97, 0.04, f"{SHAPE_V3[hab]}\nratio {r.r_val:.2f} ({r.r_lo:.2f}\u2013{r.r_hi:.2f})",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.7', alpha=0.9))
    ax.set_xlim(0, 80)
for ax in axes[-1]: ax.set_xlabel('|Latitude| (degrees)')
for ax in axes[:, 0]: ax.set_ylabel('Fitted richness')
fig.tight_layout()
fig.savefig(f'{FIG}/fig2_ldg_curves.png'); fig.savefig(f'{FIG}/fig2_ldg_curves.svg')
plt.close(fig)
print('fig2 done')

# ---------------- Fig 3: robustness dumbbells, grouped order ----------------
m = p02[['habitat','ratio_10v50']].merge(
    p03[['habitat','ratio1']], on='habitat').merge(
    p13[['habitat','ratio_with_unmapped']], on='habitat')
m = m.set_index('habitat').reindex(ORDER).reset_index()

fig, (axa, axb) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
y = np.arange(len(m))
for i, r in m.iterrows():
    axa.plot([min(r.ratio_10v50, r.ratio1), max(r.ratio_10v50, r.ratio1)], [i, i],
             color='0.6', lw=1.2, zorder=1)
axa.scatter(m.ratio_10v50, y, color='black', s=28, zorder=2, label='raw')
axa.scatter(m.ratio1, y, color='#0072B2', s=28, zorder=2, label='study random effect')
axa.axvline(1, ls='--', color='0.4', lw=0.8)
axa.set_xscale('log'); axa.set_xlim(0.15, 8)
axa.set_xticks([0.25, 0.5, 1, 2, 4]); axa.set_xticklabels(['0.25','0.5','1','2','4'])
axa.set_yticks(y); axa.set_yticklabels([HABITAT_LABELS[h] for h in m.habitat], fontsize=8.5)
axa.invert_yaxis()
axa.set_xlabel('Gradient ratio (10\u00b0 / 50\u00b0)')
axa.set_title('a  Study-composition adjustment', loc='left', fontsize=10)
axa.legend(fontsize=8, loc='lower left', frameon=False)

for i, r in m.iterrows():
    axb.plot([min(r.ratio_10v50, r.ratio_with_unmapped), max(r.ratio_10v50, r.ratio_with_unmapped)],
             [i, i], color='0.6', lw=1.2, zorder=1)
axb.scatter(m.ratio_10v50, y, color='black', s=28, zorder=2, label='raw')
axb.scatter(m.ratio_with_unmapped, y, color='#D55E00', s=28, zorder=2, label='unmapped-fraction covariate')
axb.axvline(1, ls='--', color='0.4', lw=0.8)
axb.set_xscale('log'); axb.set_xlim(0.15, 8)
axb.set_xticks([0.25, 0.5, 1, 2, 4]); axb.set_xticklabels(['0.25','0.5','1','2','4'])
axb.set_xlabel('Gradient ratio (10\u00b0 / 50\u00b0)')
axb.set_title('b  Reference-catalog bias adjustment', loc='left', fontsize=10)
axb.legend(fontsize=8, loc='lower left', frameon=False)
fig.tight_layout()
fig.savefig(f'{FIG}/fig3_robustness.png'); fig.savefig(f'{FIG}/fig3_robustness.svg')
plt.close(fig)
print('fig3 done')
