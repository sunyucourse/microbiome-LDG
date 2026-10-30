"""Coverage-based diversity standardization (Chao & Jost 2012) for MicrobeAtlas OTU tables.

Implements the same mathematics as iNEXT's coverage-based rarefaction/extrapolation:
- interpolation: exact hypergeometric expectation via gammaln
- extrapolation: Chao1-based unseen-species estimator
Reference: Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation:
standardizing samples by completeness rather than size. Ecology 93, 2533-2547.
"""
import numpy as np
from numba import njit
from math import lgamma as gammaln


@njit(cache=True)
def coverage_at_size(counts, n, m):
    """Good's coverage of a subsample of size m drawn from a sample with counts (sum=n).

    C(m) = 1 - sum_i (x_i/n) * P(species i absent in subsample of size m)
    with P(absent) the hypergeometric probability computed via gammaln.
    """
    if m >= n:
        return 1.0
    deficit = 0.0
    for x in counts:
        if x > n - m:
            p_absent = 0.0
        else:
            lg = (gammaln(n - x + 1) + gammaln(n - m + 1)
                  - gammaln(n - x - m + 1) - gammaln(n + 1))
            p_absent = np.exp(lg) if lg < 0.0 else 0.0
        deficit += (x / n) * p_absent
    return 1.0 - deficit


@njit(cache=True)
def richness_at_size(counts, n, m):
    """Expected observed richness in a subsample of size m (interpolation)."""
    if m >= n:
        return float(len(counts))
    s = 0.0
    for x in counts:
        if x > n - m:
            p_absent = 0.0
        else:
            lg = (gammaln(n - x + 1) + gammaln(n - m + 1)
                  - gammaln(n - x - m + 1) - gammaln(n + 1))
            p_absent = np.exp(lg) if lg < 0.0 else 0.0
        s += 1.0 - p_absent
    return s


@njit(cache=True)
def solve_m_for_coverage(counts, n, target_cov):
    """Find subsample size m such that coverage_at_size(m) = target_cov (bisection).
    Only valid when target_cov <= coverage at full n. Returns m (float)."""
    lo, hi = 1.0, float(n)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if coverage_at_size(counts, n, mid) < target_cov:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@njit(cache=True)
def chao1_f0(counts, n):
    """Estimated number of unobserved species (Chao1); returns (f0, f1, f2)."""
    f1 = 0
    f2 = 0
    for x in counts:
        if x == 1:
            f1 += 1
        elif x == 2:
            f2 += 1
    if f2 > 0:
        return ((n - 1.0) / n) * (f1 * f1) / (2.0 * f2), f1, f2
    else:
        return ((n - 1.0) / n) * (f1 * (f1 - 1.0)) / (2.0 * (f2 + 1.0)), f1, f2


@njit(cache=True)
def coverage_estimate(counts, n):
    """Chao-Jost adjusted Good's coverage estimate of the full sample; returns (C, f1, f2)."""
    f1 = 0
    f2 = 0
    for x in counts:
        if x == 1:
            f1 += 1
        elif x == 2:
            f2 += 1
    if f1 == 0:
        return 1.0, f1, f2
    denom = ((n - 1.0) * f1) + 2.0 * f2
    a = ((n - 1.0) * f1) / denom if denom > 0 else 0.0
    return 1.0 - (f1 / n) * a, f1, f2


@njit(cache=True)
def extrapolated_richness(counts, n, target_cov):
    """Richness extrapolated to target coverage > observed coverage (Chao & Jost 2012)."""
    f0, f1, f2 = chao1_f0(counts, n)
    if f1 == 0 or f0 <= 0:
        return float(len(counts))
    a = f1 / (n * f0 + f1)
    base = 1.0 - a
    if base <= 0.0:
        return len(counts) + f0
    # coverage at n+m: C = 1 - (f1/n) * (1-a)^m  -> solve m for target_cov
    m_star = np.log((1.0 - target_cov) * n / f1) / np.log(base)
    if m_star <= 0:
        return float(len(counts))
    return len(counts) + f0 * (1.0 - base ** m_star)


@njit(cache=True)
def standardized_richness(counts, target_cov):
    """Coverage-standardized q=0 richness at target coverage (inter- or extrapolated).

    Returns (richness_at_target, n_reads, coverage_observed)."""
    n = 0
    for x in counts:
        n += x
    if n == 0:
        return 0.0, 0, 0.0
    cov, f1, f2 = coverage_estimate(counts, n)
    if target_cov <= cov:
        m = solve_m_for_coverage(counts, n, target_cov)
        return richness_at_size(counts, n, m), n, cov
    else:
        return extrapolated_richness(counts, n, target_cov), n, cov


@njit(cache=True)
def hill_q1(counts):
    """Observed Hill number q=1 (exp Shannon)."""
    n = 0
    for x in counts:
        n += x
    if n == 0:
        return 0.0
    h = 0.0
    for x in counts:
        p = x / n
        h -= p * np.log(p)
    return np.exp(h)
