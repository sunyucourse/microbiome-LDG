#!/usr/bin/env python3
"""Figure S4: taxonomic drivers of the natural-soil reverse LDG (P17)."""
import sys
sys.path.insert(0, '/workspace/ldg')
from figstyle import *
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from scipy import stats

d = pd.read_csv('/workspace/ldg/p17_soil_phylum_drivers.csv')
gr = pd.read_csv('/workspace/ldg/p17_soil_grid_phylum_richness.csv', index_col=0)

SHORT = {'Bacteria;Verrucomicrobia':'Verrucomicrobia','Bacteria;Acidobacteria':'Acidobacteriota',
         'Bacteria;Planctomycetes':'Planctomycetes','Bacteria;Actinobacteria':'Actinobacteria',
         'Bacteria;Proteobacteria':'Proteobacteria','Bacteria;Bacteroidetes':'Bacteroidota',
         'Bacteria;Firmicutes':'Firmicutes','Bacteria;Chloroflexi':'Chloroflexi',
         'Bacteria;Gemmatimonadetes':'Gemmatimonadota','Bacteria;Cyanobacteria':'Cyanobacteria',
         'Bacteria;Nitrospirae':'Nitrospirota','Archaea;Thaumarchaeota':'Thaumarchaeota',
         'Bacteria;Armatimonadetes':'Armatimonadota','Bacteria;Elusimicrobia':'Elusimicrobia',
         'Bacteria;Fibrobacteres':'Fibrobacterota','Bacteria;Fusobacteria':'Fusobacteria',
         'Bacteria;Deinococcus-Thermus':'Deinococcus-Thermus'}
def short(p):
    if p in SHORT: return SHORT[p]
    return p.split(';')[-1] if ';' in p else p

# major phyla: present in >=300 grids, ranked by rho
maj = d[d.n_grids_present >= 300].sort_values('spearman_rho', ascending=True).reset_index(drop=True)

fig, (axa, axb) = plt.subplots(1, 2, figsize=(10, 5.2), gridspec_kw={'width_ratios':[1.15,1]})

# (a) rho per major phylum
y = np.arange(len(maj))
colors = ['#D55E00' if r > 0 else '#0072B2' for r in maj.spearman_rho]
axa.barh(y, maj.spearman_rho, color=colors, height=0.65)
axa.set_yticks(y)
axa.set_yticklabels([short(p) for p in maj.phylum], fontsize=8)
axa.axvline(0, color='black', lw=0.7)
axa.set_xlabel('Spearman rho (phylum OTU richness vs |latitude|)')
axa.set_title('a  Phylum-level drivers (grids >= 300)', loc='left', fontsize=10)
for i, (r, p) in enumerate(zip(maj.spearman_rho, maj.P)):
    star = '***' if p < 1e-4 else ('**' if p < 1e-3 else ('*' if p < 0.05 else 'n.s.'))
    x = r + (0.012 if r > 0 else -0.012)
    axa.text(x, i, star, va='center', ha='left' if r > 0 else 'right', fontsize=7)
axa.set_xlim(-0.08, 0.48)

# (b) top driver vs latitude: Acidobacteriota + Verrucomicrobia grid medians
top2 = ['Bacteria;Acidobacteria', 'Bacteria;Verrucomicrobia']
for p, colr in zip(top2, ['#009E73', '#D55E00']):
    c = gr[gr[p] > 0]
    axb.scatter(c.alat, c[p], s=9, alpha=0.35, color=colr, rasterized=True)
    rho, pv = stats.spearmanr(gr.alat, gr[p])
    xs = np.linspace(0, 80, 50)
    b = np.polyfit(gr.alat, np.log1p(gr[p]), 1)
    axb.plot(xs, np.expm1(np.polyval(b, xs)), color=colr, lw=1.6,
             label=f'{short(p)} (rho={rho:.2f})')
axb.set_xlabel('|Latitude| (degrees)')
axb.set_ylabel('OTU richness per grid (median sample)')
axb.set_yscale('symlog', linthresh=1)
axb.legend(fontsize=8, frameon=False, loc='upper left')
axb.set_title('b  Top two drivers vs latitude', loc='left', fontsize=10)

fig.tight_layout()
fig.savefig('/workspace/figures_ms/figS4_soil_drivers.png')
fig.savefig('/workspace/figures_ms/figS4_soil_drivers.svg')
print('S4 done')
