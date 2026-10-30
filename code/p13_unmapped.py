#!/usr/bin/env python3
"""P1-3: unmapped-read fraction vs latitude, per habitat, + covariate sensitivity.
Input: p13_sample_reads.tsv (MAP_SID total mapped n_otu) from refilt headers.
unmapped_fraction = 1 - mapped/total.
"""
import pandas as pd, numpy as np
from scipy import stats as st

SHARED = '/mnt/shared-workspace/shared/ldg_data'

reads = pd.read_csv('/workspace/data/processed/p13_sample_reads.tsv', sep=' ',
                    names=['MAP_SID','total','mapped','n_otu'])
reads['unmapped_frac'] = 1 - reads.mapped / reads.total
print('samples with read totals:', len(reads))
print('unmapped_frac quantiles:', reads.unmapped_frac.quantile([0,.25,.5,.75,1]).round(3).to_dict())

meta = pd.read_csv(f'{SHARED}/qc_sample_metadata.csv',
                   usecols=['MAP_SID','lat','lon','habitat_groups','bl_tier2'])
m = meta.merge(reads, on='MAP_SID')
m = m[m.lat.notna() & ~m.bl_tier2]
rows = []
for r in m.itertuples():
    if not isinstance(r.habitat_groups, str): continue
    for h in r.habitat_groups.split('|'):
        if h: rows.append((h, abs(r.lat), r.unmapped_frac))
e = pd.DataFrame(rows, columns=['habitat','lat_abs','unmapped_frac'])

# 1) per-habitat Spearman of unmapped fraction vs latitude
res = []
for hab, g in e.groupby('habitat'):
    if len(g) < 100: continue
    r = st.spearmanr(g.lat_abs, g.unmapped_frac)
    res.append(dict(habitat=hab, n=len(g), rho=round(r.statistic,3), p=round(r.pvalue,4),
                    median_frac=round(g.unmapped_frac.median(),3)))
res = pd.DataFrame(res).sort_values('rho')
print(res.to_string(index=False))
res.to_csv('/workspace/data/processed/p13_unmapped_vs_lat.csv', index=False)

# 2) grid-level covariate sensitivity: does adding unmapped fraction change gradients?
m['glat'] = np.floor(m.lat).astype(int); m['glon'] = np.floor(m.lon).astype(int)
m['grid_id'] = (m.glat + 90) * 360 + (m.glon + 180)
g = m.groupby(['habitat_groups','grid_id']).agg(lat_abs=('lat', lambda x: abs(np.floor(x.median()))),
    unmapped=('unmapped_frac','median'), n=('MAP_SID','size')).reset_index()
g = g[g.n >= 3]
g.to_csv('/workspace/data/processed/p13_grid_unmapped.csv', index=False)
print('grid-level table saved:', len(g))
