"""
module2_filtration.py — Glomerular Filtration Module

Computes GFR from Starling forces across the glomerular capillary wall.

Physics (Starling equation for ultrafiltration):
  GFR = Kf × (P_gc − P_bc − π_gc)

Where:
  P_gc = glomerular capillary hydrostatic pressure (from Module 1)
  P_bc = Bowman's capsule hydrostatic pressure (~18 mmHg, relatively constant)
  π_gc = plasma oncotic pressure in glomerular capillary (increases along length)

Oncotic pressure rises along the capillary as protein-free ultrafiltrate is
removed, concentrating plasma proteins. We use an average oncotic pressure.

References:
  - Guyton & Hall, Ch. 27: Glomerular Filtration
  - Deen et al., Am J Physiol (1972) — Starling model of glomerular filtration
  - StatPearls: Glomerular Filtration Rate physiology
"""

import numpy as np


def compute_oncotic_pressure(RPF, GFR_estimate, params):
    """
    Compute average glomerular oncotic pressure accounting for
    protein concentration along the capillary.

    As filtrate is removed, protein concentrates:
      C_protein_efferent = C_protein_afferent / (1 − FF)
    Average π ≈ mean of afferent and efferent oncotic pressure.

    π = a1 × C + a2 × C²  (Landis-Pappenheimer approximation)
    Simplified: π ≈ pi_coeff × C  (linear for moderate concentrations)
    """
    f = params.filt

    # Filtration fraction estimate
    if RPF > 0:
        FF = min(GFR_estimate / RPF, 0.45)  # physiological cap
    else:
        FF = f.FF_normal

    C_afferent = f.plasma_protein  # g/dL
    C_efferent = C_afferent / max(1.0 - FF, 0.3)

    pi_afferent = f.pi_coeff * C_afferent
    pi_efferent = f.pi_coeff * C_efferent

    # Average oncotic pressure along capillary
    pi_avg = (pi_afferent + pi_efferent) / 2.0

    return pi_avg


def compute_GFR(P_gc, RPF, params, GFR_prev=None):
    """
    Compute glomerular filtration rate.

    GFR = Kf × (P_gc − P_bc − π_gc_avg)

    Uses iterative approach since π depends on GFR (through FF),
    and GFR depends on π. We iterate until convergence.

    Parameters
    ----------
    P_gc : float — Glomerular capillary pressure (mmHg)
    RPF : float — Renal plasma flow (mL/min)
    params : Parameters
    GFR_prev : float or None — Previous GFR for initial guess

    Returns
    -------
    dict with keys:
        GFR — Glomerular filtration rate (mL/min)
        P_net — Net filtration pressure (mmHg)
        pi_gc — Average oncotic pressure (mmHg)
        FF — Filtration fraction
        filtered_Na — Filtered sodium load (mEq/min)
        filtered_water — Filtered water volume (mL/min)
    """
    f = params.filt
    tNa = params.tubNa

    # Initial GFR guess
    GFR = GFR_prev if GFR_prev is not None else f.GFR_normal

    # Iterative solution (GFR and π are coupled)
    for _ in range(10):
        pi_gc = compute_oncotic_pressure(RPF, GFR, params)
        P_net = P_gc - f.P_bc - pi_gc
        GFR_new = f.Kf * max(P_net, 0.0)
        if abs(GFR_new - GFR) < 0.01:
            break
        GFR = 0.7 * GFR_new + 0.3 * GFR  # damped update for stability

    GFR = max(GFR, 0.0)

    # Filtration fraction
    FF = GFR / RPF if RPF > 0 else 0.0

    # Filtered loads
    # Na: GFR (mL/min) × plasma [Na] (mEq/L) / 1000 (mL→L) = mEq/min
    filtered_Na = GFR * tNa.plasma_Na / 1000.0  # mEq/min

    # Water: GFR itself is the filtered water volume in mL/min
    filtered_water = GFR  # mL/min

    return {
        "GFR": GFR,
        "P_net": max(P_net, 0.0),
        "pi_gc": pi_gc,
        "FF": FF,
        "filtered_Na": filtered_Na,
        "filtered_water": filtered_water,
    }
