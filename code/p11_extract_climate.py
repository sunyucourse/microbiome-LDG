#!/usr/bin/env python3
"""Extract per-grid climate: WorldClim BIO1/4/12 (land) + WOA23 SST (ocean)."""
import rasterio, numpy as np, pandas as pd, xarray as xr

gt = pd.read_csv('/mnt/results/tables/grid_table.csv')
cells = gt[['habitat','grid_id','glat','glon']].drop_duplicates()
print('cells:', len(cells))

# ---- WorldClim: average 6x6 pixels within each 1-deg cell ----
def wc_extract(path):
    src = rasterio.open(path)
    out = {}
    for glat, glon in zip(cells.glat, cells.glon):
        win = rasterio.windows.from_bounds(glon, glat, glon+1, glat+1, src.transform)
        win = win.round_offsets().round_lengths()
        try:
            a = src.read(1, window=win, masked=True)
            v = float(a.mean()) if a.count() else np.nan
        except Exception:
            v = np.nan
        out[(glat, glon)] = v
    src.close()
    return out

bio1 = wc_extract('/workspace/wc/wc2.1_10m_bio_1.tif')    # degC x10
bio4 = wc_extract('/workspace/wc/wc2.1_10m_bio_4.tif')    # SD x100
bio12 = wc_extract('/workspace/wc/wc2.1_10m_bio_12.tif')  # mm
cells['bio1'] = [bio1[(a,b)] for a,b in zip(cells.glat, cells.glon)]        # degC (verified empirically)
cells['bio4'] = [bio4[(a,b)]/100 for a,b in zip(cells.glat, cells.glon)]    # SD x100 -> degC
cells['bio12'] = [bio12[(a,b)] for a,b in zip(cells.glat, cells.glon)]      # mm

# ---- WOA23 SST (0 m), 1-deg; nearest valid cell within 2 deg for coastal NaN ----
woa = xr.open_dataset('/workspace/woa23_sst.nc', decode_times=False)
sst = woa['t_an'].isel(time=0, depth=0).values  # (180, 360)
lats = woa.lat.values; lons = woa.lon.values
def sst_at(glat, glon):
    i = int(np.argmin(np.abs(lats - (glat + 0.5))))
    j = int(np.argmin(np.abs(lons - (glon + 0.5))))
    v = sst[i, j]
    if np.isfinite(v): return float(v)
    for r in (1, 2):
        sub = sst[max(0,i-r):i+r+1, max(0,j-r):j+r+1]
        if np.isfinite(sub).any():
            return float(np.nanmean(sub))
    return np.nan
cells['sst'] = [sst_at(a, b) for a, b in zip(cells.glat, cells.glon)]

MARINE = {'marine_planktonic','marine_sediment'}
cells['temp'] = np.where(cells.habitat.isin(MARINE), cells.sst, cells.bio1)
cells.to_csv('/workspace/data/processed/p11_grid_climate.csv', index=False)
print(cells.groupby('habitat')[['bio1','sst','temp']].apply(lambda g: g.temp.notna().sum()).to_string())
print(cells[['bio1','bio4','bio12','sst','temp']].describe().round(2))
