"""Odds ratios, exact CIs and Mantel-Haenszel stratified ORs (no statsmodels needed)."""
import math

from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio


def or_ci(a, b, c, d):
    """2x2: [[a exposed-case, b exposed-noncase],[c unexposed-case, d unexposed-noncase]].
    Returns (sample OR, conditional-MLE OR, 95% CI low, high, two-sided Fisher p)."""
    table = [[a, b], [c, d]]
    _, p = fisher_exact(table)
    res = odds_ratio(table, kind="conditional")
    lo, hi = res.confidence_interval(confidence_level=0.95)
    sample = float("inf") if (b * c) == 0 else (a * d) / (b * c)
    return sample, res.statistic, lo, hi, p


def woolf_ci(a, b, c, d):
    """Haldane-Anscombe corrected log-OR CI (for cells that may be 0)."""
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    lor = math.log((a_ * d_) / (b_ * c_))
    se = math.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    return math.exp(lor), math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se), se


def mantel_haenszel(strata):
    """strata: iterable of (a, b, c, d) per stratum.
    Returns (OR_MH, ci_lo, ci_hi, chi2_MH, p, n_strata_used, n_total)."""
    from scipy.stats import chi2 as chi2dist
    num = den = 0.0
    # Robins-Breslow-Greenland variance pieces
    P_R = P_S = Q_R = Q_S = R = S = 0.0
    used = 0
    tot = 0
    # MH chi-square pieces
    o_sum = e_sum = v_sum = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        tot += n
        if (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        used += 1
        Ri = a * d / n
        Si = b * c / n
        num += Ri
        den += Si
        R += Ri
        S += Si
        Pi = (a + d) / n
        Qi = (b + c) / n
        P_R += Pi * Ri
        P_S += Pi * Si
        Q_R += Qi * Ri
        Q_S += Qi * Si
        o_sum += a
        e_sum += (a + b) * (a + c) / n
        if n > 1:
            v_sum += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    if den == 0 or num == 0:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), used, tot
    or_mh = num / den
    var = P_R / (2 * R * R) + (P_S + Q_R) / (2 * R * S) + Q_S / (2 * S * S)
    se = math.sqrt(var)
    lo = math.exp(math.log(or_mh) - 1.96 * se)
    hi = math.exp(math.log(or_mh) + 1.96 * se)
    if v_sum > 0:
        chi2 = (abs(o_sum - e_sum) - 0.5) ** 2 / v_sum
        p = chi2dist.sf(chi2, 1)
    else:
        chi2, p = float("nan"), float("nan")
    return or_mh, lo, hi, chi2, p, used, tot


def wilson(k, n):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    z = 1.959964
    p = k / n
    c = z * z / n
    centre = (p + c / 2) / (1 + c)
    half = z * math.sqrt(p * (1 - p) / n + c / (4 * n)) / (1 + c)
    return p, max(0.0, centre - half), min(1.0, centre + half)
