#!/usr/bin/env python3
"""H3: engineered/urban habitats decoupled from latitude, replaced by
population/GDP gradients.

Steps:
1. Take wastewater + urban_built grid tables (from main grid_table.csv),
   attach country (modal iso3 of constituent samples) and World Bank
   indicators (population, GDP per capita, sanitation).
2. Variance decomposition of log(rich90): |lat| vs log pop vs log GDPpc
   vs sanitation — sequential + partial R2 via OLS (statsmodels).
3. Same-latitude contrast: for each engineered grid, natural soil/freshwater
   grids within +-2 deg latitude — compare diversity residuals.

Inputs:
  grid_table.csv (main, per-habitat grid rows with grid_id)
  qc_sample_metadata.csv (MAP_SID -> iso3, habitat_groups)
  sample_diversity.csv
  worldbank_recent.csv (iso3, pop, gdp_pc, sanit_basic, sanit_safe)
Outputs:
  h3_engineered_grids.csv, h3_variance_decomposition.csv, h3_contrast.csv
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

SHARED = '/mnt/shared-workspace/shared/ldg_data'
OUT = '/workspace/data/processed'
WB = '/workspace/data/processed/worldbank_recent.csv'
os.makedirs(OUT, exist_ok=True)

ENGINEERED = ['wastewater', 'urban_built']
NATURAL_REF = ['soil_natural', 'freshwater']

def main():
    gt = pd.read_csv(f'{OUT}/grid_table.csv')  # per-habitat alpha grid table
    print('grid_table cols:', gt.columns.tolist())

    # modal iso3 per grid from sample metadata
    meta = pd.read_csv(f'{SHARED}/qc_sample_metadata.csv',
                       usecols=['MAP_SID', 'lat', 'lon', 'iso3', 'bl_tier2'])
    meta = meta[(~meta['bl_tier2'].fillna(True)) & meta['lat'].notna() & meta['iso3'].notna()]
    meta['grid'] = (np.clip(np.floor(meta['lat']).astype(int) + 90, 0, 179) * 360
                    + np.clip(np.floor(meta['lon']).astype(int) + 180, 0, 359))
    modal_iso = (meta.groupby('grid')['iso3']
                 .agg(lambda x: x.mode().iloc[0] if len(x.mode()) else np.nan))

    wb = pd.read_csv(WB)
    wb = wb.rename(columns={'iso3c': 'iso3'})

    eng = gt[gt['habitat'].isin(ENGINEERED)].copy()
    eng['iso3'] = eng['grid_id'].map(modal_iso)
    eng = eng.merge(wb, on='iso3', how='left')
    eng.to_csv(f'{OUT}/h3_engineered_grids.csv', index=False)
    print(f'engineered grids: {len(eng)}, with WB data: {eng["gdp_pc"].notna().sum()}')

    # ---- variance decomposition ----
    results = []
    for hab in ENGINEERED + NATURAL_REF:
        sub = gt[gt['habitat'] == hab].copy()
        if hab in ENGINEERED:
            sub = eng[eng['habitat'] == hab].copy()
            sub['log_pop'] = np.log10(sub['pop'])
            sub['log_gdp'] = np.log10(sub['gdp_pc'])
            sub['sanit'] = sub['sanit_basic']
            predictors = ['lat_abs', 'log_pop', 'log_gdp', 'sanit']
        else:
            predictors = ['lat_abs']
        sub = sub.dropna(subset=['richness', 'lat_abs'] + [p for p in predictors if p != 'lat_abs'])
        sub = sub[np.isfinite(sub['richness']) & (sub['richness'] > 0)]
        if len(sub) < 30:
            print(f'{hab}: too few grids ({len(sub)})'); continue
        y = np.log(sub['richness'])
        # full model
        X = sm.add_constant(sub[predictors].astype(float))
        full = sm.OLS(y, X).fit()
        row = {'habitat': hab, 'n_grids': len(sub), 'R2_full': full.rsquared}
        # leave-one-out partial R2
        for p in predictors:
            others = [q for q in predictors if q != p]
            if others:
                Xr = sm.add_constant(sub[others].astype(float))
                r2_red = sm.OLS(y, Xr).fit().rsquared
            else:
                r2_red = 0.0
            row[f'partial_R2_{p}'] = max(0.0, full.rsquared - r2_red)
            row[f'beta_{p}'] = full.params.get(p, np.nan)
            row[f'p_{p}'] = full.pvalues.get(p, np.nan)
        results.append(row)
    vd = pd.DataFrame(results)
    vd.to_csv(f'{OUT}/h3_variance_decomposition.csv', index=False)
    print(vd.to_string())

    # ---- same-latitude contrast (glat is already signed: grid centre = glat + 0.5) ----
    gt['lat_signed'] = gt['glat'] + 0.5
    nat = gt[gt['habitat'].isin(NATURAL_REF)].copy()
    rows = []
    for hab in ENGINEERED:
        e = gt[gt['habitat'] == hab]
        for _, r in e.iterrows():
            m = nat[(nat['lat_signed'] >= r['lat_signed'] - 2) & (nat['lat_signed'] <= r['lat_signed'] + 2)]
            if len(m) >= 3:
                rows.append({'habitat': hab, 'grid_id': r['grid_id'],
                             'lat': r['lat_signed'], 'rich_eng': r['richness'],
                             'rich_nat_median': m['richness'].median(),
                             'n_nat_grids': len(m)})
    con = pd.DataFrame(rows)
    if len(con):
        con['log_ratio'] = np.log(con['rich_eng'] / con['rich_nat_median'])
        con.to_csv(f'{OUT}/h3_contrast.csv', index=False)
        print(con.groupby('habitat')['log_ratio'].describe().to_string())

if __name__ == '__main__':
    main()
