#!/usr/bin/env python3
"""P17: taxonomic drivers of the natural-soil reverse LDG.
Per-sample per-phylum OTU richness (mapped, non-organelle reads; coverage>=0.90 samples),
grid medians (>=3 samples/1-deg cell), Spearman + effort-partial correlations vs |latitude|.
Run on machine 'taxon' (32GB)."""
import h5py, numpy as np, pandas as pd
from scipy import stats

BIOM = '/workspace/data/companion.biom'
TAX = '/workspace/data/taxonomic_mapping.tsv'
SOIL = '/mnt/shared-workspace/shared/soil_natural_samples_v4.tsv'
OUT = '/mnt/shared-workspace/shared/'

# ---------- load BIOM structure ----------
f = h5py.File(BIOM, 'r')
sids = np.array([x.decode() for x in f['sample/ids'][:]])
oids = np.array([x.decode() for x in f['observation/ids'][:]])
print('samples', len(sids), 'otus', len(oids))
print('obs matrix keys:', list(f['observation/matrix'].keys()))
print('sample matrix keys:', list(f['sample/matrix'].keys()) if 'sample/matrix' in f else 'none')

# ---------- taxonomy mapping ----------
tax = pd.read_csv(TAX, sep='\t', dtype=str)
print('tax columns:', list(tax.columns), len(tax))
oid_index = pd.Index(oids)
best_col, best_hit = None, 0
for col in tax.columns:
    if col == 'Taxonomy':
        continue
    hit = oid_index.intersection(pd.Index(tax[col].dropna().unique())).size
    print(col, 'hits:', hit)
    if hit > best_hit:
        best_col, best_hit = col, hit
print('using id column:', best_col)
tmap = tax.drop_duplicates(best_col).set_index(best_col).Taxonomy
o2p_full = pd.Series(oids).map(tmap)
# phylum = first two taxonomy levels (domain;phylum)
o2p = o2p_full.fillna('').str.split(';').str[:2].agg(';'.join).where(o2p_full.notna())
print('mapped to taxonomy:', o2p.notna().sum(), 'of', len(oids))

is_unmapped = pd.Series(oids).str.contains('nmapped', case=False).values
is_organelle = pd.Series(oids).str.startswith('E97').values | \
               o2p_full.fillna('').str.contains('hloroplast|itochondri', case=False).values
print('unmapped rows:', is_unmapped.sum(), 'organelle:', is_organelle.sum())
keep_otu = ~(is_unmapped | is_organelle) & o2p.notna().values

ph_names = np.sort(o2p[keep_otu].unique())
ph2code = {p: i for i, p in enumerate(ph_names)}
otu_phcode = np.full(len(oids), -1, dtype=np.int32)
otu_phcode[np.where(keep_otu)[0]] = o2p[keep_otu].map(ph2code).values

# ---------- soil samples present in BIOM ----------
soil = pd.read_csv(SOIL, sep='\t')
sid2col = {s: i for i, s in enumerate(sids)}
soil['col'] = soil.MAP_SID.map(sid2col)
soil = soil.dropna(subset=['col'])
soil['col'] = soil.col.astype(int)
print('soil samples in BIOM:', len(soil))

# ---------- read matrix (CSR by observation) ----------
data = f['observation/matrix/data'][:]
indices = f['observation/matrix/indices'][:]
indptr = f['observation/matrix/indptr'][:]
print('nnz:', len(data))
wanted = np.zeros(len(sids), dtype=bool)
wanted[soil.col.values] = True
rows_otu = np.repeat(np.arange(len(oids)), np.diff(indptr))
mask = wanted[indices]
sel_otu, sel_sam, sel_cnt = rows_otu[mask], indices[mask], data[mask]
print('selected nnz:', len(sel_cnt))
order = np.argsort(sel_sam, kind='stable')
sel_sam, sel_otu, sel_cnt = sel_sam[order], sel_otu[order], sel_cnt[order]

# ---------- per-sample coverage + per-phylum richness ----------
n_ph = len(ph_names)
res = np.zeros((len(soil), n_ph), dtype=np.int32)
cov_ok = np.zeros(len(soil), dtype=bool)
nreads = np.zeros(len(soil), dtype=np.int64)
cols = soil.col.values
for j in range(len(soil)):
    c = cols[j]
    i0 = np.searchsorted(sel_sam, c); i1 = np.searchsorted(sel_sam, c, side='right')
    if i1 == i0:
        continue
    otus_c = sel_otu[i0:i1]; cnts_c = sel_cnt[i0:i1]
    mapped_mask = ~is_unmapped[otus_c]
    n = int(cnts_c[mapped_mask].sum())
    if n == 0:
        continue
    f1 = int((cnts_c[mapped_mask] == 1).sum())
    if 1 - f1 / n < 0.90:
        continue
    cov_ok[j] = True; nreads[j] = n
    keep = otu_phcode[otus_c] >= 0
    res[j] = np.bincount(otu_phcode[otus_c[keep]], minlength=n_ph)
print('coverage-pass samples:', cov_ok.sum())

soil['cov_ok'] = cov_ok; soil['n_reads'] = nreads
soil = soil[soil.cov_ok].reset_index(drop=True)
res = res[cov_ok]
print('after coverage filter:', len(soil))

# ---------- grid aggregation (>=3 samples) ----------
res_df = pd.DataFrame(res, columns=ph_names)
res_df['grid_id'] = soil.grid_id.values
grid_n = soil.groupby('grid_id').size()
ok = grid_n[grid_n >= 3].index
res_df = res_df[res_df.grid_id.isin(ok)]
gridmed = res_df.groupby('grid_id').median()
grid_n = grid_n.loc[gridmed.index]
grid_alat = (soil.drop_duplicates('grid_id').set_index('grid_id').glat + 0.5).abs().loc[gridmed.index]
print('grids >=3 samples:', len(gridmed))

# ---------- per-phylum correlations ----------
out = []
logn = np.log(grid_n.values)
alat = grid_alat.values
for p in ph_names:
    v = gridmed[p].values
    if (v > 0).sum() < 30:
        continue
    rho, pv = stats.spearmanr(alat, v)
    rx = stats.rankdata(alat); ry = stats.rankdata(v); rc = stats.rankdata(logn)
    bx = np.polyfit(rc, rx, 1); by = np.polyfit(rc, ry, 1)
    pr, pp = stats.pearsonr(rx - np.polyval(bx, rc), ry - np.polyval(by, rc))
    lo = alat <= 23.5; hi = alat >= 50
    trop = np.median(v[lo]) if lo.sum() > 5 else np.nan
    high = np.median(v[hi]) if hi.sum() > 5 else np.nan
    out.append([p, int((v > 0).sum()), rho, pv, pr, pp, trop, high,
                (high / trop) if trop and trop > 0 else np.nan, np.median(v)])
res_tab = pd.DataFrame(out, columns=['phylum', 'n_grids_present', 'spearman_rho', 'P',
                                     'partial_r_effort', 'P_partial', 'median_tropical',
                                     'median_highlat', 'ratio_high_over_trop', 'median_per_grid'])
res_tab = res_tab.sort_values('spearman_rho', ascending=False)
res_tab.to_csv(OUT + 'p17_soil_phylum_drivers.csv', index=False)
gridmed.assign(alat=alat, n=grid_n.values).to_csv(OUT + 'p17_soil_grid_phylum_richness.csv')
pd.set_option('display.width', 250)
print(res_tab.head(15).to_string(index=False))
print('...bottom:')
print(res_tab.tail(8).to_string(index=False))
