#!/usr/bin/env python3
"""H1 meta-regression: does LDG strength track climate exposure across habitats?
Reads fit_ldg_gam.R strength output, joins habitat exposure scores, computes
Spearman correlation + leave-one-habitat-out cross-validation of the sign/rank.
Outputs h1_meta_regression.csv and prints summary.
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, '/workspace/ldg')
from figstyle import EXPOSURE, HABITAT_LABELS

STRENGTH = sys.argv[1] if len(sys.argv) > 1 else '/workspace/data/processed/gam_all_strength.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else '/workspace/data/processed/h1_meta_regression.csv'

s = pd.read_csv(STRENGTH)
s['exposure'] = s['habitat'].map(EXPOSURE)
s = s.dropna(subset=['exposure'])
from scipy.stats import t as tdist
s['lat2_p'] = 2 * tdist.sf(np.abs(s['lat2_t']), df=s['n_grids'] - 5)
print(s[['habitat', 'n_grids', 'lat2_coef', 'lat2_p', 'eq_vs_highlat_ratio', 'exposure']].to_string())

# primary metric: equator/high-latitude ratio (log scale); secondary: lat2 coefficient
for metric in ['eq_vs_highlat_ratio', 'lat2_coef']:
    x = s['exposure'].to_numpy(float)
    y = s[metric].to_numpy(float)
    if metric == 'eq_vs_highlat_ratio':
        y = np.log(y)
    rho, p = spearmanr(x, y)
    print(f'\n{metric}: Spearman rho = {rho:.3f}, p = {p:.4f}')

    # leave-one-habitat-out CV: predict sign of relationship from n-1 habitats
    signs = []
    for i in range(len(s)):
        m = np.ones(len(s), bool); m[i] = False
        r_loo, _ = spearmanr(x[m], y[m])
        signs.append(np.sign(r_loo) == np.sign(rho))
    print(f'  LOO sign stability: {np.mean(signs):.2%}')

s.to_csv(OUT, index=False)
print('saved', OUT)
