"""
module1_hemodynamics.py — Renal Hemodynamics Module

Computes renal blood flow (RBF), renal plasma flow (RPF), and glomerular
capillary pressure (P_gc) from arterial pressure, arteriolar resistances,
and modulatory signals (RSNA, Ang II, TGF).

Physics:
  RBF = (MAP − P_venous) / RVR
  RVR = R_aa + R_ea
  P_gc = MAP − RBF × R_aa

Autoregulation is achieved through:
  1. Myogenic response (stretch-dependent afferent arteriolar tone)
  2. Tubuloglomerular feedback (TGF signal modulates R_aa)
  3. RSNA modulation of both arterioles
  4. Ang II preferentially constricts efferent arteriole

References:
  - Guyton & Hall, Ch. 26–27
  - Karaaslan et al. (2005), Eqs. for renal vascular resistance
  - Hallow et al. (2014), renal hemodynamic submodel
"""

import numpy as np


def compute_afferent_resistance(MAP, RSNA, angII, TGF_signal, params):
    """
    Compute afferent arteriolar resistance.

    R_aa = R_aa_base × (1 + myogenic) × (1 + TGF) × (1 + RSNA_effect) × (1 + AngII_effect)

    Parameters
    ----------
    MAP : float — Mean arterial pressure (mmHg)
    RSNA : float — Normalized renal sympathetic nerve activity (1.0 = normal)
    angII : float — Normalized angiotensin II (1.0 = normal)
    TGF_signal : float — TGF-mediated fractional change in R_aa
    params : Parameters
    """
    h = params.hemo
    hr = params.horm

    R_aa_base = h.RVR_normal * h.R_aa_fraction  # baseline R_aa

    # Myogenic response: constriction when MAP rises, dilation when MAP falls
    # Sigmoidal within autoregulatory range
    MAP_deviation = MAP - h.MAP_normal
    myogenic_factor = h.myogenic_gain * MAP_deviation
    myogenic_factor = np.clip(myogenic_factor, -0.3, 0.5)

    # TGF: comes from module3/TGF computation
    tgf_factor = TGF_signal  # already computed as fractional change

    # RSNA: sympathetic vasoconstriction
    rsna_factor = hr.RSNA_afferent_sensitivity * (RSNA - 1.0)

    # Ang II: mild afferent constriction at high levels
    angII_factor = hr.angII_afferent_sensitivity * max(angII - 1.5, 0.0)

    # Combined multiplicative modulation
    R_aa = R_aa_base * (1.0 + myogenic_factor) * (1.0 + tgf_factor) * \
           (1.0 + rsna_factor) * (1.0 + angII_factor)

    return max(R_aa, 0.01)  # prevent non-physical values


def compute_efferent_resistance(RSNA, angII, params):
    """
    Compute efferent arteriolar resistance.

    Ang II preferentially constricts the efferent arteriole (Guyton).
    RSNA also contributes to efferent tone.

    R_ea = R_ea_base × (1 + AngII_effect) × (1 + RSNA_effect)
    """
    h = params.hemo
    hr = params.horm

    R_ea_base = h.RVR_normal * h.R_ea_fraction

    # Ang II: strong efferent constriction (key RAAS mechanism to maintain GFR)
    angII_factor = hr.angII_efferent_sensitivity * (angII - 1.0)

    # RSNA: milder efferent effect
    rsna_factor = hr.RSNA_efferent_sensitivity * (RSNA - 1.0)

    R_ea = R_ea_base * (1.0 + angII_factor) * (1.0 + rsna_factor)

    return max(R_ea, 0.01)


def compute_renal_hemodynamics(MAP, RSNA, angII, TGF_signal, params):
    """
    Main hemodynamics computation.

    Returns
    -------
    dict with keys:
        RBF — Renal blood flow (mL/min)
        RPF — Renal plasma flow (mL/min)
        RVR — Renal vascular resistance (mmHg·min/mL)
        R_aa — Afferent resistance
        R_ea — Efferent resistance
        P_gc — Glomerular capillary pressure (mmHg)
    """
    h = params.hemo

    R_aa = compute_afferent_resistance(MAP, RSNA, angII, TGF_signal, params)
    R_ea = compute_efferent_resistance(RSNA, angII, params)

    RVR = R_aa + R_ea
    P_venous = 4.0  # Renal venous pressure (mmHg) — approximately constant

    # Ohm's law for flow
    RBF = max((MAP - P_venous) / RVR, 10.0)  # floor to prevent negative

    # Plasma flow
    RPF = RBF * (1.0 - h.Hct)

    # Glomerular capillary pressure
    # P_gc = MAP - RBF × R_aa  (pressure drop across afferent arteriole)
    P_gc = MAP - RBF * R_aa
    P_gc = max(P_gc, 20.0)  # physiological floor

    return {
        "RBF": RBF,
        "RPF": RPF,
        "RVR": RVR,
        "R_aa": R_aa,
        "R_ea": R_ea,
        "P_gc": P_gc,
    }
