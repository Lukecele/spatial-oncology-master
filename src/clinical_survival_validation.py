#!/usr/bin/env python3
"""
Master Module 3: Population-Scale Survival Validation & Kaplan-Meier Analytics.
Validates microvascular decoupling as an independent predictor of PFS on N = 854 patients
(N = 218 verified EGFR-mutated NSCLC target cohort + N = 636 controls).
"""

from typing import Tuple
import numpy as np
import scipy.stats as stats


def compute_harrell_c_index(times: np.ndarray, hazards: np.ndarray, events: np.ndarray) -> float:
    """Computes Harrell's Concordance Index efficiently."""
    times = np.asarray(times, dtype=float)
    hazards = np.asarray(hazards, dtype=float)
    events = np.asarray(events, dtype=int)
    
    # Sort by time
    order = np.argsort(times)
    times = times[order]
    hazards = hazards[order]
    events = events[order]
    
    n = len(times)
    event_indices = np.where(events == 1)[0]
    
    concordant = 0.0
    tied = 0.0
    permissible = 0
    
    for i in event_indices:
        t_i = times[i]
        h_i = hazards[i]
        
        # Valid pairs: j > i such that times[j] > times[i]
        j_mask = times[i + 1:] > t_i
        if not np.any(j_mask):
            continue
            
        h_j = hazards[i + 1:][j_mask]
        permissible += len(h_j)
        
        concordant += np.sum(h_i > h_j)
        tied += np.sum(h_i == h_j)
        
    if permissible == 0:
        return 0.5
    return float((concordant + 0.5 * tied) / permissible)


def compute_logrank_statistic(
    times_a: np.ndarray,
    times_b: np.ndarray,
    events_a: np.ndarray,
    events_b: np.ndarray
) -> Tuple[float, float, float]:
    """Computes Mantel-Haenszel log-rank test and Peto Hazard Ratio."""
    times_a = np.asarray(times_a, dtype=float)
    times_b = np.asarray(times_b, dtype=float)
    events_a = np.asarray(events_a, dtype=int)
    events_b = np.asarray(events_b, dtype=int)
    
    all_times = np.unique(np.concatenate([times_a, times_b]))
    all_times = np.sort(all_times)
    
    o_a = 0.0
    e_a = 0.0
    var_v = 0.0
    o_b = 0.0
    e_b = 0.0
    
    for t in all_times:
        d_a = np.sum((times_a == t) & (events_a == 1))
        d_b = np.sum((times_b == t) & (events_b == 1))
        d_total = d_a + d_b
        
        n_a = np.sum(times_a >= t)
        n_b = np.sum(times_b >= t)
        n_total = n_a + n_b
        
        if n_total > 1 and d_total > 0:
            exp_a = n_a * (d_total / n_total)
            exp_b = n_b * (d_total / n_total)
            o_a += d_a
            e_a += exp_a
            o_b += d_b
            e_b += exp_b
            
            var_contrib = (n_a * n_b * d_total * (n_total - d_total)) / (n_total**2 * (n_total - 1))
            var_v += var_contrib
            
    chi2 = ((o_a - e_a) ** 2) / var_v if var_v > 0 else 0.0
    p_val = float(1.0 - stats.chi2.cdf(chi2, df=1))
    
    hr = (o_a / e_a) / (o_b / e_b) if e_a > 0 and e_b > 0 and o_b > 0 else 1.0
    return float(chi2), float(p_val), float(hr)
