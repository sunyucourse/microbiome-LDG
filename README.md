# Habitat-specific latitudinal diversity gradients in the global microbiome

Analysis code and derived data for the manuscript:

> Sun Y, Fu X. Habitat-specific latitudinal diversity gradients in the global microbiome. (submitted)

## Overview

We reanalyzed 1,039,362 quality-controlled microbiome samples from the MicrobeAtlas
database (352,941 samples eligible after coverage, geographic, and habitat filters),
spanning twelve habitat types, and estimated habitat-specific latitudinal diversity
gradients (LDGs) under a single uniform pipeline. Four direction types emerged
(classic decline, reverse increase, mid-latitude peak, flat); environmental mediation
was habitat-specific, consistent with a limiting-factor covariance hypothesis.

## Data source

Raw data are **not** redistributed here. Download from MicrobeAtlas:

- OTU table (coverage-filtered companion BIOM, 1,039,362 samples):
  https://microbeatlas.org/downloads/paper/metag_covfilter/samples-otus.97.metag.minfilter.minCov90.noMulticell.rod2025companion.biom.gz
- Sample metadata and taxonomy mapping:
  https://microbeatlas.org/downloads/paper/taxonomic_mapping.tsv
- Database website: https://microbeatlas.org/

Environmental layers used: WorldClim 2.1 (BIO1, BIO4, BIO12), World Ocean Atlas 2023
(sea-surface temperature, nitrate), SoilGrids v2.0 (soil pH, soil organic carbon),
GHS-POP (population), World Bank country indicators (GDP per capita, sanitation).

## Repository structure

```
code/       analysis scripts (Python + R)
data/       derived data tables (grid-level and summary results)
figures/    final manuscript figures (PNG + SVG)
```

## Analysis pipeline (run order)

1. `code/diversity.py` — coverage-standardized richness (rich90; Chao & Jost 2012)
   from the BIOM table; sample eligibility filters.
2. `code/taxonomy.py` — OTU taxonomy mapping; exclusion of unmapped and organelle OTUs.
3. `code/aggregate_grids.py` — habitat assignment (two-tier keyword blacklists) and
   aggregation to 1-degree grid cells (>= 3 samples per cell).
4. `code/fit_ldg_gam.R` — per-habitat GAMs of log richness vs absolute latitude
   (mgcv, REML); 10/50-degree ratios with posterior-simulation CIs.
5. `code/p02_gradient_metrics.R` — direction typology and curve-shape classification.
6. `code/p03_study_re.R` — study random-effect (study-composition) adjustment.
7. `code/p13_unmapped.py`, `code/p13_unmapped_sens.R` — reference-catalog bias
   (unmapped-fraction) sensitivity.
8. `code/p11_extract_climate.py`, `code/p11_temperature.R` — climate-layer extraction
   and temperature sensitivity.
9. `code/p14_soil_ph.py`, `code/p14_soil_soc.py` — soil pH / organic-carbon mediation.
10. `code/p16_tost.R` — host thermoregulation equivalence test (TOST).
11. `code/p17_soil_taxon_drivers.py` — phylum-level drivers of the soil reverse gradient.
12. `code/sensitivity.py` — robustness suite (rich95, Hill q=1, 2-degree grid,
    blacklist tiers, hemisphere, technology).
13. `code/h1_meta_regression.py`, `code/h2_analysis.py`, `code/h3_analysis.py` —
    climate-exposure hypothesis and supporting analyses.
14. `code/fig*.py`, `code/figs_v3.py`, `code/supp_figs_v3.py` — figure generation.

## Software versions

- R 4.4.2; mgcv 1.9.4
- Python 3.11; pandas 2.3.1, scipy 1.15.0, statsmodels 0.14.6, numpy, matplotlib, seaborn

## License

Code is released under the MIT License (see LICENSE). Derived data tables are released
under CC-BY 4.0.

## Citation

<!-- TODO: add Zenodo DOI badge after archiving, e.g.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->
