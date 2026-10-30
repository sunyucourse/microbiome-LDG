#!/usr/bin/env python3
"""H2: endotherm vs ectotherm host LDG.
Builds grid-level diversity tables per (body_site x thermo) group from
sample-level diversity + curated host metadata.

Inputs (shared workspace):
  sample_diversity.csv     MAP_SID + diversity stats (Unmapped/E97 excluded)
  animal_host_metadata.csv MAP_SID, is_gut/is_skin, thermo, host_class, sub_envs,
                           lat/lon/Technology/bl_tier2
Outputs (OUT dir):
  grid_h2_<site>_<thermo>[_wildonly].csv   grid tables for GAM
"""
import os
import numpy as np
import pandas as pd

SHARED = '/mnt/shared-workspace/shared/ldg_data'
OUT = '/workspace/data/processed'
os.makedirs(OUT, exist_ok=True)
MIN_SAMPLES = 3

def build_grid_table(df, group_name):
    """Emit grid table in the schema fit_ldg_gam.R expects."""
    df = df.copy()
    df['glat'] = np.clip(np.floor(df['lat']).astype(int) + 90, 0, 179)
    df['glon'] = np.clip(np.floor(df['lon']).astype(int) + 180, 0, 359)
    df['grid_id'] = df['glat'] * 360 + df['glon']
    g = df.groupby('grid_id')
    out = g.agg(glat=('glat', 'first'), glon=('glon', 'first'),
                n_samples=('rich90', 'size'),
                lat_abs=('lat', lambda x: x.abs().median()),
                richness=('rich90', 'median'), rich95=('rich95', 'median'),
                hill_q1=('hill_q1', 'median'),
                s_obs=('s_obs', 'median'),
                tech_frac_wgs=('Technology', lambda x: (x == 'WGS').mean()))
    out = out[out['n_samples'] >= MIN_SAMPLES].reset_index()
    out['habitat'] = group_name
    out['depth_m'] = np.nan
    return out[['habitat', 'grid_id', 'glat', 'glon', 'lat_abs', 'richness',
                'rich95', 'hill_q1', 's_obs', 'n_samples', 'depth_m', 'tech_frac_wgs']]

def main():
    sdiv = pd.read_csv(f'{SHARED}/sample_diversity.csv')
    sdiv = sdiv[sdiv['coverage'] >= 0.90]
    host = pd.read_csv(f'{SHARED}/animal_host_metadata.csv',
                       usecols=['MAP_SID', 'is_gut', 'is_skin', 'thermo', 'host_class',
                                'sub_envs', 'lat', 'lon', 'Technology', 'bl_tier2'])
    df = host.merge(sdiv, on='MAP_SID')
    df = df[(~df['bl_tier2'].fillna(True)) & df['lat'].notna() & df['thermo'].notna()]
    df['site'] = np.where(df['is_gut'], 'gut', np.where(df['is_skin'], 'skin', ''))
    df = df[df['site'] != '']
    print(f'H2 samples with coords+diversity: {len(df):,}')
    print(df.groupby(['site', 'thermo']).size())

    for site in ['gut', 'skin']:
        parts = []
        for thermo in ['endotherm', 'ectotherm']:
            sub = df[(df['site'] == site) & (df['thermo'] == thermo)]
            gt = build_grid_table(sub, f'{site}_{thermo}')
            parts.append(gt)
            print(f'{site}/{thermo}: {len(gt)} grids, {len(sub)} samples')
        fn = f'{OUT}/grid_h2_{site}.csv'
        pd.concat(parts).to_csv(fn, index=False)
        print(f'{fn} written')

    # sensitivity: exclude human-commensal / lab / livestock hosts by sub_env keywords
    pat = (r'mouse|murine|\brat\b|pig|swine|cow|cattle|bovine|chicken|dog|canine|'
           r'cat|feline|sheep|ovine|goat|horse|equine|macaque|human')
    df2 = df[~df['sub_envs'].fillna('').str.lower().str.contains(pat, regex=True)]
    print(f'wild-only samples: {len(df2):,}')
    for site in ['gut', 'skin']:
        parts = []
        for thermo in ['endotherm', 'ectotherm']:
            sub = df2[(df2['site'] == site) & (df2['thermo'] == thermo)]
            parts.append(build_grid_table(sub, f'{site}_{thermo}'))
        fn = f'{OUT}/grid_h2_{site}_wildonly.csv'
        pd.concat(parts).to_csv(fn, index=False)
        print(f'{fn} written')

if __name__ == '__main__':
    main()
