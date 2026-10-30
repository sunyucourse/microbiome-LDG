#!/usr/bin/env python3
"""Aggregate per-sample diversity + COO shards into a grid-level analysis table.

Inputs (from /workspace/out/chunk_*.npz on heavy-1):
  - per-sample stats (BIOM sample order)
  - COO triplets (pair_id, otu_idx, count) for analysis samples
Output: grid_table.csv with per (habitat, grid) alpha-median and pooled-gamma diversity.
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, '/workspace/ldg')
from diversity import standardized_richness

OUT = '/workspace/out'
HABITATS = ['soil_natural', 'soil_agricultural', 'marine_planktonic', 'marine_sediment',
            'freshwater', 'wastewater', 'urban_built', 'plant_rhizosphere',
            'plant_phyllosphere', 'animal_gut', 'animal_skin', 'human_gut']
NGRID = 180 * 360


def load_chunks():
    stats_list, sid_list, cp, co, cc = [], [], [], [], []
    for fp in sorted(glob.glob(f'{OUT}/chunk_*.npz')):
        z = np.load(fp, allow_pickle=False)
        stats_list.append(z['stats'])
        sid_list.append(z['sids'])
        cp.append(z['coo_pair']); co.append(z['coo_otu']); cc.append(z['coo_cnt'])
    stats = np.concatenate(stats_list)
    sids = np.concatenate(sid_list)
    return stats, sids, np.concatenate(cp), np.concatenate(co), np.concatenate(cc)


def main():
    print('loading chunks...', flush=True)
    stats, sids, coo_pair, coo_otu, coo_cnt = load_chunks()
    cols = ['n_reads', 's_obs', 'f1', 'f2', 'coverage', 'rich90', 'rich95', 'hill_q1']
    sdf = pd.DataFrame(stats, columns=cols)
    sdf['MAP_SID'] = sids
    sdf = sdf[sdf['n_reads'] > 0]
    # paper's minCov90 criterion, now computable on mapped reads (Unmapped excluded):
    sdf = sdf[sdf['coverage'] >= 0.90]
    print('sample stats:', len(sdf), flush=True)
    sdf.to_csv(f'{OUT}/sample_diversity.csv', index=False)

    # metadata for grid assignment
    meta = pd.read_csv('/workspace/data/qc_sample_metadata.csv',
                       dtype={'lat': 'float32', 'lon': 'float32'})
    meta = meta.set_index('MAP_SID')
    sdf = sdf.join(meta[['lat', 'lon', 'Technology', 'bl_tier1', 'bl_tier2', 'habitat_groups']],
                   on='MAP_SID')
    # depth covariate (sample-level parsed measurements, keyed by sample accession)
    try:
        meas = pd.read_csv('/workspace/data/samples.measurements.tsv', sep='\t', dtype=str)
        dep = meas[meas['Meta_variable'] == 'depth_m'][['SampleID', 'Value']].copy()
        dep['depth_m'] = pd.to_numeric(dep['Value'], errors='coerce')
        dep = dep.dropna(subset=['depth_m']).groupby('SampleID')['depth_m'].median()
        sdf['SampleAcc'] = sdf['MAP_SID'].str.split('.').str[-1]
        sdf['depth_m'] = sdf['SampleAcc'].map(dep)
    except Exception as e:
        print('depth merge failed:', e, flush=True)
        sdf['depth_m'] = np.nan
    sdf = sdf[sdf['lat'].notna() & sdf['habitat_groups'].fillna('').ne('')]
    sdf['glat'] = np.floor(sdf['lat']).astype(int)
    sdf['glon'] = np.floor(sdf['lon']).astype(int)
    sdf['grid_id'] = (sdf['glat'] + 90) * 360 + (sdf['glon'] + 180)
    sdf['lat_abs'] = sdf['lat'].abs()
    sdf['is_wgs'] = (sdf['Technology'] == 'WGS').astype(float)

    # ---- alpha (sample-median) grid table, tier2-filtered ----
    rows = []
    sdf['bl_tier2'] = sdf['bl_tier2'].fillna(False).astype(bool)
    sdf['bl_tier1'] = sdf['bl_tier1'].fillna(False).astype(bool)
    sdf_t2 = sdf[~sdf['bl_tier2']]
    exp = sdf_t2.assign(habitat=sdf_t2['habitat_groups'].str.split(r'\|')).explode('habitat')
    for (hab, gid), g in exp.groupby(['habitat', 'grid_id']):
        if len(g) < 3:
            continue
        rows.append(dict(habitat=hab, grid_id=int(gid),
                         glat=int(g['glat'].iloc[0]), glon=int(g['glon'].iloc[0]),
                         lat_abs=float(g['lat_abs'].median()),
                         richness=float(g['rich90'].median()),
                         rich95=float(g['rich95'].median()),
                         hill_q1=float(g['hill_q1'].median()),
                         s_obs=float(g['s_obs'].median()),
                         n_samples=len(g),
                         depth_m=float(g['depth_m'].median()) if g['depth_m'].notna().any() else np.nan,
                         tech_frac_wgs=float(g['is_wgs'].mean())))
    alpha_tab = pd.DataFrame(rows)
    print('alpha grid table:', alpha_tab.shape, flush=True)

    # ---- pooled gamma per (habitat, grid) from COO ----
    if len(coo_pair) > 0:
        pair_max = int(coo_pair.max()) + 1
        notu = 133249
        print(f'building pair x OTU matrix: {len(coo_pair)} nnz, {pair_max} pairs', flush=True)
        M = sparse.csr_matrix((coo_cnt, (coo_pair, coo_otu)), shape=(pair_max, notu))
        gamma_rows = []
        nz_pairs = np.unique(coo_pair)
        for pid in nz_pairs:
            hab_code = pid // NGRID
            gid = pid % NGRID
            counts = M.getrow(pid).data.astype(np.int64)
            if counts.sum() == 0:
                continue
            r90, n, cov = standardized_richness(counts, 0.90)
            r95, _, _ = standardized_richness(counts, 0.95)
            gamma_rows.append(dict(pair_id=int(pid), habitat=HABITATS[hab_code], grid_id=int(gid),
                                   gamma_rich90=float(r90), gamma_rich95=float(r95),
                                   gamma_reads=int(n), gamma_cov=float(cov),
                                   gamma_sobs=int(len(counts))))
        gamma_tab = pd.DataFrame(gamma_rows)
        print('gamma grid table:', gamma_tab.shape, flush=True)
        grid_tab = alpha_tab.merge(gamma_tab.drop(columns=['pair_id']), on=['habitat', 'grid_id'], how='left')
    else:
        grid_tab = alpha_tab
    grid_tab.to_csv(f'{OUT}/grid_table.csv', index=False)
    print('GRID_TABLE_DONE', flush=True)


if __name__ == '__main__':
    main()
