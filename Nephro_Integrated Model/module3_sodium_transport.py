"""
module3_sodium_transport.py — Tubular Sodium Transport Module

Computes segmental sodium reabsorption and excretion along:
  1. Proximal Tubule (PT): ~67% of filtered load
  2. Loop of Henle (LoH): ~25% of delivered load (thick ascending limb)
  3. Distal Tubule (DT): ~5% of delivered load
  4. Collecting Duct (CD): ~2-5% of delivered load (fine-tuning)

Modulation:
  - PT: Ang II (↑), RSNA (↑), pressure natriuresis (↓ with ↑ MAP)
  - LoH: relatively fixed (loop diuretics target here)
  - DT: aldosterone (↑), thiazides target here
  - CD: aldosterone (↑), ANP (↓)

Also computes macula densa sodium delivery for TGF feedback.

References:
  - Guyton & Hall, Ch. 28–29
  - Karaaslan et al. (2005) — sodium reabsorption equations
  - Hallow et al. (2014) — tubular transport submodel
"""

import numpy as np


def compute_PT_reabsorption(filtered_Na, MAP, angII, RSNA, params):
    """
    Proximal tubule sodium reabsorption.

    Modulated by:
      - Glomerulotubular balance (reabsorption scales with filtered load — built in)
      - Ang II stimulation
      - RSNA stimulation
      - Pressure natriuresis (↑ MAP → ↓ reabsorption via peritubular Starling forces)

    Returns: (Na_reabsorbed_PT, Na_delivered_to_LoH)
    """
    t = params.tubNa

    # Baseline fractional reabsorption
    frac = t.frac_reabs_PT

    # Ang II modulation: increases reabsorption
    angII_effect = t.PT_angII_sensitivity * (angII - 1.0)

    # RSNA modulation: increases reabsorption
    rsna_effect = t.PT_rsna_sensitivity * (RSNA - 1.0)

    # Pressure natriuresis: ↑ MAP → ↓ reabsorption
    # This is a key long-term mechanism (Guyton)
    pressure_effect = -t.PT_pressure_sensitivity * (MAP - params.hemo.MAP_normal)

    # Net fractional reabsorption (clamped to physiological range)
    frac_modulated = frac + angII_effect + rsna_effect + pressure_effect
    frac_modulated = np.clip(frac_modulated, 0.40, 0.85)

    Na_reabsorbed = filtered_Na * frac_modulated
    Na_delivered_to_LoH = filtered_Na - Na_reabsorbed

    return Na_reabsorbed, Na_delivered_to_LoH, frac_modulated


def compute_LoH_reabsorption(Na_delivered, params, loop_diuretic_factor=1.0):
    """
    Loop of Henle sodium reabsorption (thick ascending limb).

    Relatively fixed fraction; target of loop diuretics.
    loop_diuretic_factor: 1.0 = normal, <1.0 = diuretic effect

    Returns: (Na_reabsorbed_LoH, Na_delivered_to_DT)
    """
    t = params.tubNa

    frac = t.frac_reabs_LoH * loop_diuretic_factor
    frac = np.clip(frac, 0.0, 0.95)

    Na_reabsorbed = Na_delivered * frac
    Na_delivered_to_DT = Na_delivered - Na_reabsorbed

    return Na_reabsorbed, Na_delivered_to_DT


def compute_DT_reabsorption(Na_delivered, aldosterone, params, thiazide_factor=1.0):
    """
    Distal tubule sodium reabsorption.

    Modulated by aldosterone. Target of thiazide diuretics.

    Returns: (Na_reabsorbed_DT, Na_delivered_to_CD)
    """
    t = params.tubNa

    frac = t.frac_reabs_DT
    aldo_effect = t.DT_aldosterone_sensitivity * (aldosterone - 1.0)

    frac_modulated = (frac + aldo_effect) * thiazide_factor
    frac_modulated = np.clip(frac_modulated, 0.0, 0.80)

    Na_reabsorbed = Na_delivered * frac_modulated
    Na_delivered_to_CD = Na_delivered - Na_reabsorbed

    return Na_reabsorbed, Na_delivered_to_CD


def compute_CD_reabsorption(Na_delivered, aldosterone, ANP, params,
                            K_sparing_factor=1.0):
    """
    Collecting duct sodium reabsorption.

    Fine-tuning segment. Modulated by:
      - Aldosterone (↑ reabsorption)
      - ANP (↓ reabsorption, promotes natriuresis)
      - K-sparing diuretics reduce reabsorption here

    Returns: (Na_reabsorbed_CD, Na_excreted)
    """
    t = params.tubNa

    frac = t.frac_reabs_CD
    aldo_effect = t.CD_aldosterone_sensitivity * (aldosterone - 1.0)
    anp_effect = -t.CD_anp_sensitivity * (ANP - 1.0)

    frac_modulated = (frac + aldo_effect + anp_effect) * K_sparing_factor
    frac_modulated = np.clip(frac_modulated, 0.0, 0.98)

    Na_reabsorbed = Na_delivered * frac_modulated
    Na_excreted = Na_delivered - Na_reabsorbed

    return Na_reabsorbed, Na_excreted


def compute_TGF_signal(Na_delivered_to_DT, params):
    """
    Tubuloglomerular feedback signal.

    Macula densa senses NaCl delivery at the early distal tubule.
    When delivery increases → TGF constricts afferent arteriole → reduces GFR.
    When delivery decreases → TGF dilates afferent arteriole → increases GFR.

    Signal modeled as sigmoidal function of normalized macula densa Na delivery.

    Returns: TGF_signal — fractional change in afferent resistance
             Positive = constriction, Negative = dilation
    """
    tgf = params.tgf

    # Normal macula densa delivery (mEq/min)
    # At baseline: filtered_Na × (1 - frac_PT) × (1 - frac_LoH)
    t = params.tubNa
    normal_delivery = (t.filtered_Na_normal / 1440.0) * (1 - t.frac_reabs_PT) * \
                      (1 - t.frac_reabs_LoH)

    # Normalized delivery
    if normal_delivery > 0:
        x = Na_delivered_to_DT / normal_delivery
    else:
        x = 1.0

    # Sigmoidal TGF response
    # When x > 1 (high delivery): constriction (positive signal)
    # When x < 1 (low delivery): dilation (negative signal)
    sigmoid = 1.0 / (1.0 + np.exp(-tgf.TGF_slope * (x - tgf.TGF_midpoint)))

    # Map sigmoid (0–1) to range (−max_dilation, +max_constriction)
    TGF_range = tgf.TGF_max_constriction + tgf.TGF_max_dilation
    TGF_signal = -tgf.TGF_max_dilation + sigmoid * TGF_range

    return TGF_signal


def compute_sodium_transport(filtered_Na, MAP, angII, RSNA, aldosterone, ANP, params,
                             drug_effects=None):
    """
    Full tubular sodium transport computation.

    Parameters
    ----------
    filtered_Na : float — Filtered sodium load (mEq/min)
    MAP : float — Mean arterial pressure (mmHg)
    angII : float — Normalized angiotensin II
    RSNA : float — Normalized RSNA
    aldosterone : float — Normalized aldosterone
    ANP : float — Normalized ANP
    params : Parameters
    drug_effects : dict or None — Drug modulation factors

    Returns
    -------
    dict with segmental reabsorption and excretion values
    """
    if drug_effects is None:
        drug_effects = {}

    loop_factor = drug_effects.get("loop_diuretic", 1.0)
    thiazide_factor = drug_effects.get("thiazide", 1.0)
    K_sparing_factor = drug_effects.get("K_sparing", 1.0)

    # Segment-by-segment computation
    Na_reabs_PT, Na_to_LoH, frac_PT = compute_PT_reabsorption(
        filtered_Na, MAP, angII, RSNA, params)

    Na_reabs_LoH, Na_to_DT = compute_LoH_reabsorption(
        Na_to_LoH, params, loop_factor)

    Na_reabs_DT, Na_to_CD = compute_DT_reabsorption(
        Na_to_DT, aldosterone, params, thiazide_factor)

    Na_reabs_CD, Na_excreted = compute_CD_reabsorption(
        Na_to_CD, aldosterone, ANP, params, K_sparing_factor)

    # TGF signal based on macula densa delivery (= Na_to_DT)
    TGF_signal = compute_TGF_signal(Na_to_DT, params)

    # Total reabsorption
    Na_total_reabs = Na_reabs_PT + Na_reabs_LoH + Na_reabs_DT + Na_reabs_CD
    frac_total = Na_total_reabs / filtered_Na if filtered_Na > 0 else 0.0

    return {
        "Na_reabs_PT": Na_reabs_PT,
        "Na_reabs_LoH": Na_reabs_LoH,
        "Na_reabs_DT": Na_reabs_DT,
        "Na_reabs_CD": Na_reabs_CD,
        "Na_to_LoH": Na_to_LoH,
        "Na_to_DT": Na_to_DT,          # Macula densa delivery
        "Na_to_CD": Na_to_CD,
        "Na_excreted": Na_excreted,     # Final urinary Na excretion (mEq/min)
        "Na_total_reabs": Na_total_reabs,
        "frac_total_reabs": frac_total,
        "frac_PT": frac_PT,
        "TGF_signal": TGF_signal,
    }
