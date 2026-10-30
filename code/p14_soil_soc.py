#!/usr/bin/env python3
"""Query SoilGrids pH (0-5cm) for soil/rhizosphere grid cells; save grid pH table."""
import pandas as pd, numpy as np, urllib.request, json, time, sys

gt = pd.read_csv('/mnt/results/tables/grid_table.csv')
cells = gt[gt.habitat.isin(['soil_natural','soil_agricultural','plant_rhizosphere','plant_phyllosphere'])][
    ['habitat','grid_id','glat','glon']].drop_duplicates()
print('cells to query:', len(cells))

def query_ph(lat, lon, tries=3):
    url = (f'https://rest.isric.org/soilgrids/v2.0/properties/query?lon={lon}&lat={lat}'
           f'&property=soc&depth=0-5cm&value=mean')
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                j = json.load(r)
            v = j['properties']['layers'][0]['depths'][0]['values']['mean']
            return v / 10 if v is not None else np.nan
        except Exception:
            time.sleep(1.5)
    return np.nan

out = []
for i, r in enumerate(cells.itertuples()):
    out.append(query_ph(r.glat + 0.5, r.glon + 0.5))
    if (i + 1) % 100 == 0:
        print(f'{i+1}/{len(cells)}', flush=True)
    time.sleep(0.12)
cells['soil_soc'] = out
cells.to_csv('/workspace/data/processed/p14_soil_soc.csv', index=False)
print('done, non-null:', cells.soil_soc.notna().sum())
