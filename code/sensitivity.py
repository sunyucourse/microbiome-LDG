#!/usr/bin/env python3
"""P8 sensitivity suite: re-aggregate alpha-diversity grid tables under
alternative analysis choices and run the same GAM per habitat.

Scenarios:
  base        1-deg grids, >=3 samples, tier2 blacklist, rich90   (reference)
  grid2deg    2-deg grids, >=5 samples
  nh_only     Northern Hemisphere grids only
  hill_q1     Hill q=1 as response instead of coverage-standardized richness
  tier1       tier1 blacklist instead of tier2
  amplicon    AMPLICON samples only
  wgs         WGS samples only

Outputs: /workspace/data/processed/sens_<scenario>.csv grid tables
         (GAMs run separately via fit_ldg_gam.R per scenario)
"""
import os, sys
import numpy as np
import pandas as pd

SHARED = '/mnt/shared-workspace/shared/ldg_data'
OUT = '/workspace/data/processed'
os.makedirs(OUT, exist_ok=True)

def aggregate(df, value_col, deg=1, min_n=3):
    df = df.copy()
    df['glat'] = np.clip(np.floor(df['lat'] / deg).astype(int) * deg + 90, 0, 179)
    df['glon'] = np.clip(np.floor(df['lon'] / deg).astype(int) * deg + 180, 0, 359)
    df['grid_id'] = (df['glat'] // deg) * (360 // deg) + (df['glon'] // deg)
    g = df.groupby(['habitat', 'grid_id'])
    out = g.agg(glat=('glat', 'first'), glon=('glon', 'first'),
                n_samples=(value_col, 'size'),
                lat_abs=('lat', lambda x: x.abs().median()),
                richness=(value_col, 'median'),
                s_obs=('s_obs', 'median'),
                depth_m=('depth_m', 'median'),
                tech_frac_wgs=('Technology', lambda x: (x == 'WGS').mean()))
    out = out[out['n_samples'] >= min_n].reset_index()
    out['depth_m'] = out['depth_m'].where(out.groupby('habitat')['depth_m']
                                          .transform(lambda s: s.notna().sum() >= 30))
    return out

def main():
    sdiv = pd.read_csv(f'{SHARED}/sample_diversity.csv')
    sdiv = sdiv[sdiv['coverage'] >= 0.90]
    meta = pd.read_csv(f'{SHARED}/qc_sample_metadata.csv',
                       usecols=['MAP_SID', 'lat', 'lon', 'Technology', 'bl_tier1',
                                'bl_tier2', 'habitat_groups'])
    df = sdiv.merge(meta, on='MAP_SID')
    df = df[df['lat'].notna() & df['habitat_groups'].notna()]
    # depth from measurements (long format: SampleID, Meta_variable, Value)
    try:
        meas = pd.read_csv(f'{SHARED}/samples.measurements.tsv', sep='\t', dtype=str)
        dep = meas[meas['Meta_variable'] == 'depth_m'][['SampleID', 'Value']].copy()
        dep['depth_m'] = pd.to_numeric(dep['Value'], errors='coerce')
        dep = dep.dropna(subset=['depth_m']).groupby('SampleID')['depth_m'].median()
        df['SampleAcc'] = df['MAP_SID'].str.split('.').str[-1]
        df['depth_m'] = df['SampleAcc'].map(dep)
    except Exception as e:
        print('depth merge failed:', e)
        df['depth_m'] = np.nan
    df = df.drop(columns=['SampleAcc'], errors='ignore')
    df['habitat_groups'] = df['habitat_groups'].str.split(r'\|')
    df = df.explode('habitat_groups').rename(columns={'habitat_groups': 'habitat'})
    print(f'base pool: {len(df):,} sample-habitat rows')

    scenarios = {
        'grid2deg':  dict(deg=2, min_n=5),
        'nh_only':   dict(deg=1, min_n=3, pre=lambda d: d[d['lat'] >= 0]),
        'hill_q1':   dict(deg=1, min_n=3, value='hill_q1'),
        'tier1':     dict(deg=1, min_n=3, bl='bl_tier1'),
        'amplicon':  dict(deg=1, min_n=3, pre=lambda d: d[d['Technology'] == 'AMPLICON']),
        'wgs':       dict(deg=1, min_n=3, pre=lambda d: d[d['Technology'] == 'WGS']),
    }
    for name, kw in scenarios.items():
        d = df[~df[kw.get('bl', 'bl_tier2')].fillna(True)]
        if 'pre' in kw:
            d = kw['pre'](d)
        gt = aggregate(d, kw.get('value', 'rich90'), deg=kw['deg'], min_n=kw['min_n'])
        fn = f'{OUT}/sens_{name}.csv'
        gt.to_csv(fn, index=False)
        print(f'{name}: {len(gt)} grids -> {fn}')

if __name__ == '__main__':
    main()
