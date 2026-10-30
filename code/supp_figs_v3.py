#!/usr/bin/env python3
"""Supplementary figures for LDG manuscript v3.
S1: sensitivity of gradient ratios across variants. S2: host thermoregulation curves.
S3: soil pH / SOC vs latitude and richness (natural soils)."""
import sys
sys.path.insert(0, '/workspace/ldg')
from figstyle import *
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

OUT = '/workspace/figures_ms'
ORDER = ['soil_natural','soil_agricultural','freshwater','marine_planktonic','marine_sediment',
         'wastewater','urban_built','plant_rhizosphere','plant_phyllosphere',
         'human_gut','animal_gut','animal_skin']

# ---------------- Figure S1: sensitivity of ratios ----------------
s = pd.read_csv('/mnt/results/tables/sensitivity_summary.csv').set_index('habitat').reindex(ORDER)
variants = [('base_ratio','Base (1 deg, tier2, rich90)','black'),
            ('ratio_grid2deg','2-degree grid','#0072B2'),
            ('ratio_nh_only','Northern Hemisphere only','#56B4E9'),
            ('ratio_hill_q1','Hill q=1 (exp-Shannon)','#009E73'),
            ('ratio_tier1','Tier-1 blacklist','#E69F00'),
            ('ratio_amplicon','Amplicon only','#CC6677'),
            ('ratio_wgs','WGS only','#332288')]
fig, ax = plt.subplots(figsize=(8, 5))
y = np.arange(len(s))
for k, (col, lab, colr) in enumerate(variants):
    off = (k - 3) * 0.11
    ax.scatter(s[col], y + off, s=26, color=colr, label=lab, zorder=3)
ax.axvline(1, ls='--', color='0.4', lw=0.8)
ax.set_xscale('log'); ax.set_xlim(0.1, 12)
ax.set_xticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.125','0.25','0.5','1','2','4','8'])
ax.set_yticks(y); ax.set_yticklabels([HABITAT_LABELS[h] for h in ORDER], fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel('Gradient ratio (10\u00b0 / 50\u00b0), log scale')
ax.legend(fontsize=7.5, loc='lower left', frameon=False, ncol=2)
fig.tight_layout()
fig.savefig(f'{OUT}/figS1_sensitivity.png'); fig.savefig(f'{OUT}/figS1_sensitivity.svg')
plt.close(fig)
print('S1 done')

# ---------------- Figure S2: host thermoregulation ----------------
gut = pd.read_csv('/workspace/data/processed/gam_h2_gut_curves.csv')
skin = pd.read_csv('/workspace/data/processed/gam_h2_skin_curves.csv')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6), sharex=True)
for ax, d, ttl in ((ax1, gut, 'a  Animal gut'), (ax2, skin, 'b  Animal skin/surface')):
    for grp, colr, lab in (('endotherm','#D55E00','Endotherm hosts'),
                           ('ectotherm','#0072B2','Ectotherm hosts')):
        c = d[d.habitat.str.endswith(grp)]
        if len(c) == 0: continue
        ax.fill_between(c.lat_abs, np.exp(c.fit - 1.96*c.se), np.exp(c.fit + 1.96*c.se),
                        color=colr, alpha=0.2, linewidth=0)
        ax.plot(c.lat_abs, np.exp(c.fit), color=colr, lw=1.8, label=lab)
    ax.set_yscale('log'); ax.set_xlim(0, 75)
    ax.set_xlabel('|Latitude| (degrees)'); ax.set_title(ttl, loc='left', fontsize=10)
ax1.set_ylabel('Fitted richness')
ax1.legend(fontsize=8, frameon=False, loc='upper right')
fig.tight_layout()
fig.savefig(f'{OUT}/figS2_host_thermo.png'); fig.savefig(f'{OUT}/figS2_host_thermo.svg')
plt.close(fig)
print('S2 done')

# ---------------- Figure S3: soil pH / SOC ----------------
gt = pd.read_csv('/mnt/results/tables/grid_table.csv')
gs = gt[gt.habitat == 'soil_natural'][['grid_id','richness']].copy()
ph = pd.read_csv('/workspace/data/processed/p14_soil_ph.csv')
ph = ph[ph.habitat == 'soil_natural'][['grid_id','glat','soil_ph']]
soc = pd.read_csv('/workspace/data/processed/p14_soil_soc.csv')
soc = soc[soc.habitat == 'soil_natural'][['grid_id','soil_soc']]
m = ph.merge(soc, on='grid_id').merge(gs, on='grid_id').dropna()
m['alat'] = m.glat.abs()

fig, axes = plt.subplots(2, 2, figsize=(9, 7))
ax = axes[0,0]
ax.scatter(m.alat, m.soil_ph, s=10, alpha=0.4, color='#009E73', rasterized=True)
z = np.polyfit(m.alat, m.soil_ph, 2); xs = np.linspace(0, 75, 100)
ax.plot(xs, np.polyval(z, xs), color='black', lw=1.5)
ax.set_xlabel('|Latitude| (degrees)'); ax.set_ylabel('Soil pH (0-5 cm)')
ax.set_title('a  Soil pH vs latitude (hump-shaped)', loc='left', fontsize=10)

ax = axes[0,1]
ax.scatter(m.soil_ph, m.richness, s=10, alpha=0.4, color='#009E73', rasterized=True)
z = np.polyfit(m.soil_ph, np.log(m.richness), 2); xs = np.linspace(m.soil_ph.min(), m.soil_ph.max(), 100)
ax.plot(xs, np.exp(np.polyval(z, xs)), color='black', lw=1.5)
ax.set_yscale('log'); ax.set_xlabel('Soil pH (0-5 cm)'); ax.set_ylabel('Grid mean richness (rich90)')
ax.set_title('b  Richness vs soil pH (peak near neutral)', loc='left', fontsize=10)

ax = axes[1,0]
ax.scatter(m.alat, m.soil_soc, s=10, alpha=0.4, color='#8B5A2B', rasterized=True)
ax.set_yscale('log'); ax.set_xlabel('|Latitude| (degrees)'); ax.set_ylabel('Soil organic carbon (g/kg)')
ax.set_title('c  Soil organic carbon vs latitude', loc='left', fontsize=10)

ax = axes[1,1]
ax.scatter(m.soil_soc, m.richness, s=10, alpha=0.4, color='#8B5A2B', rasterized=True)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Soil organic carbon (g/kg)'); ax.set_ylabel('Grid mean richness (rich90)')
ax.set_title('d  Richness vs soil organic carbon (weak)', loc='left', fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUT}/figS3_soil_ph_soc.png'); fig.savefig(f'{OUT}/figS3_soil_ph_soc.svg')
plt.close(fig)
print('S3 done, n grids =', len(m))
