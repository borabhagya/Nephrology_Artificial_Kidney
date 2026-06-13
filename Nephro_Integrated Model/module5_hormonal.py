"""
module5_hormonal.py — Hormonal Regulation Module

Computes the dynamic concentrations of:
  - Renin (juxtaglomerular cells)
  - Angiotensin II (RAAS cascade)
  - Aldosterone (adrenal cortex)
  - ADH / Vasopressin (posterior pituitary)
  - ANP / Atrial Natriuretic Peptide (atrial cardiomyocytes)
  - RSNA (renal sympathetic nerve activity, normalized)

Each hormone follows first-order secretion-clearance kinetics:
  dH/dt = secretion_rate − clearance_rate × H

Secretion rates depend on physiological stimuli (MAP, Na delivery,
plasma osmolality, blood volume, etc.) forming the feedback loops.

References:
  - Karaaslan et al. (2005) — RSNA equations, renin secretion model
  - Hallow et al. (2014) — RAAS model with Ang II negative feedback on renin
  - Kutumova et al. (2021) — Calibrated agent-based BP regulation model
  - StatPearls — RAAS physiology, Vasopressin physiology
  - Bankir (2017) — ADH half-life ~3 min, osmotic threshold 280 mOsm/kg
  - StatPearls — ANP half-life 2–5 min
"""

import numpy as np


def compute_RSNA(MAP, BV, params):
    """
    Compute renal sympathetic nerve activity (normalized).

    RSNA is driven by arterial baroreceptors and cardiopulmonary receptors.
    ↓ MAP → ↑ RSNA (baroreflex)
    ↓ BV → ↑ RSNA (cardiopulmonary reflex)

    Modeled as a sigmoidal inverse relationship with MAP.

    Based on Karaaslan et al. (2005), Eqs. 2–4.
    """
    h = params.horm
    hemo = params.hemo
    vp = params.volp

    # Baroreflex component
    MAP_deviation = (MAP - hemo.MAP_normal) / hemo.MAP_normal
    baro_component = -h.RSNA_MAP_sensitivity * MAP_deviation * 100.0

    # Cardiopulmonary reflex component
    BV_deviation = (BV - hemo.BV_normal) / hemo.BV_normal
    cardiopulm_component = -h.RSNA_RAP_sensitivity * BV_deviation * 100.0

    # Total RSNA (normalized, 1.0 = normal)
    RSNA = h.RSNA_normal + baro_component + cardiopulm_component

    # Physiological bounds
    RSNA = np.clip(RSNA, 0.2, 4.0)

    return RSNA


def compute_renin_secretion_rate(MAP, Na_to_DT, RSNA, angII, ANP, params):
    """
    Compute renin secretion rate from juxtaglomerular cells.

    Stimuli for renin release (Guyton, StatPearls):
      1. ↓ Renal perfusion pressure → ↑ renin (baroreceptor mechanism)
      2. ↓ NaCl delivery to macula densa → ↑ renin
      3. ↑ RSNA (β1-adrenergic) → ↑ renin
      4. ↑ Ang II → ↓ renin (short-loop negative feedback)
      5. ↑ ANP → ↓ renin

    Returns: renin_secretion (normalized rate)
    """
    h = params.horm
    hemo = params.hemo
    t = params.tubNa

    # Normal macula densa delivery for reference
    normal_MD_delivery = (t.filtered_Na_normal / 1440.0) * (1 - t.frac_reabs_PT) * \
                         (1 - t.frac_reabs_LoH)

    # 1. Renal perfusion pressure (inverse sigmoid)
    MAP_ratio = MAP / hemo.MAP_normal
    # Inverse: high MAP suppresses, low MAP stimulates
    pressure_effect = 1.0 / (1.0 + np.exp(5.0 * (MAP_ratio - 1.0)))
    pressure_effect = 2.0 * pressure_effect  # scale so normal ≈ 1.0

    # 2. Macula densa sensing (inverse: low delivery → high renin)
    if normal_MD_delivery > 0:
        MD_ratio = Na_to_DT / normal_MD_delivery
    else:
        MD_ratio = 1.0
    MD_effect = 1.0 / (1.0 + np.exp(4.0 * (MD_ratio - 1.0)))
    MD_effect = 2.0 * MD_effect

    # 3. RSNA (direct stimulation)
    rsna_effect = 1.0 + h.renin_rsna_sensitivity * (RSNA - 1.0)
    rsna_effect = max(rsna_effect, 0.3)

    # 4. Ang II negative feedback (kicks in above normal)
    angII_inhibition = 1.0 / (1.0 + h.renin_angII_feedback * max(angII - 1.0, 0.0))

    # 5. ANP inhibition (kicks in above normal)
    anp_inhibition = 1.0 / (1.0 + h.ANP_renin_inhibition * max(ANP - 1.0, 0.0))

    # Combined secretion rate (multiplicative interactions)
    renin_secretion = h.renin_secretion_rate_normal * \
                      pressure_effect * MD_effect * rsna_effect * \
                      angII_inhibition * anp_inhibition

    return max(renin_secretion, 0.01)


def compute_renin_derivatives(renin, renin_secretion, params):
    """
    dRenin/dt = k_clear × (renin_secretion − Renin)

    At steady state with renin_secretion=1.0 and Renin=1.0, derivative=0.
    k_clearance = ln(2) / t1/2
    """
    h = params.horm
    k_clear = np.log(2) / h.renin_half_life

    d_renin = k_clear * (renin_secretion - renin)
    return d_renin


def compute_angII(renin, params, ACE_inhibition=1.0):
    """
    Angiotensin II production (quasi-steady-state due to very short half-life).

    Ang II is produced from Ang I via ACE. Since t1/2 < 1 min,
    Ang II tracks renin essentially instantaneously.

    AngII ∝ renin × ACE_activity

    ACE_inhibition: 1.0 = normal, 0.0 = full ACE inhibitor effect
    """
    h = params.horm
    angII = h.angII_from_renin_gain * renin * ACE_inhibition
    return max(angII, 0.01)


def compute_aldosterone_derivatives(aldosterone, angII, ANP, params):
    """
    dAldosterone/dt = k_clear × (secretion_target − Aldosterone)

    Secretion target determined by Ang II stimulus and ANP inhibition.
    t1/2 ≈ 20 min (hepatic clearance).
    """
    h = params.horm
    k_clear = np.log(2) / h.aldo_half_life

    # Secretion target (normalized: 1.0 at normal Ang II and ANP)
    angII_stimulus = angII  # Proportional to Ang II (1:1 at baseline)
    anp_inhibition = 1.0 / (1.0 + h.aldo_ANP_inhibition * max(ANP - 1.0, 0.0))

    target = angII_stimulus * anp_inhibition

    d_aldo = k_clear * (target - aldosterone)
    return d_aldo


def compute_ADH_derivatives(ADH, plasma_osm, BV, angII, params):
    """
    dADH/dt = k_clear × (target − ADH)

    ADH target driven by:
      1. ↑ Plasma osmolality above normal (primary)
      2. ↓ Blood volume (secondary, requires >10-20% decrease)
      3. ↑ Ang II (modest effect)

    t1/2 ≈ 4 min (Bankir 2017).
    """
    h = params.horm
    k_clear = np.log(2) / h.ADH_half_life

    # 1. Osmotic stimulus relative to MODEL's equilibrium osmolality
    # Normal osm = factor × (total_Na_normal / (ECFV_normal/1000))
    osm_normal = params.volp.Na_to_osm_factor * (
        params.volp.total_body_Na_normal / (params.volp.ECFV_normal / 1000.0))
    osm_deviation = plasma_osm - osm_normal
    osm_stimulus = h.ADH_osm_sensitivity * osm_deviation / 10.0

    # 2. Volume stimulus (requires significant hypovolemia)
    BV_ratio = BV / params.hemo.BV_normal
    if BV_ratio < 0.9:
        volume_stimulus = h.ADH_volume_sensitivity * (0.9 - BV_ratio) * 10.0
    else:
        volume_stimulus = 0.0

    # 3. Ang II stimulus
    angII_stimulus = h.ADH_angII_sensitivity * max(angII - 1.0, 0.0)

    # Target ADH (=1.0 at normal osm, normal BV, normal AngII)
    target = max(h.ADH_normal + osm_stimulus + volume_stimulus + angII_stimulus, 0.01)

    d_ADH = k_clear * (target - ADH)
    return d_ADH


def compute_ANP_derivatives(ANP, BV, params):
    """
    dANP/dt = k_clear × (target − ANP)

    ANP target set by atrial stretch (proportional to blood volume).
    t1/2 ≈ 3 min (StatPearls: 2–5 min).
    """
    h = params.horm
    k_clear = np.log(2) / h.ANP_half_life

    # Atrial stretch stimulus
    BV_ratio = BV / params.hemo.BV_normal
    target = h.ANP_normal * (1.0 + h.ANP_volume_sensitivity * (BV_ratio - 1.0))
    target = max(target, 0.1)

    d_ANP = k_clear * (target - ANP)
    return d_ANP


def compute_hormonal_derivatives(state, MAP, BV, plasma_osm, Na_to_DT, params,
                                 drug_effects=None):
    """
    Compute time derivatives for all hormonal state variables.

    Parameters
    ----------
    state : dict with current values of renin, angII, aldosterone, ADH, ANP
    MAP : float — Mean arterial pressure
    BV : float — Blood volume
    plasma_osm : float — Plasma osmolality
    Na_to_DT : float — Sodium delivery to macula densa / distal tubule
    params : Parameters
    drug_effects : dict or None

    Returns
    -------
    dict of derivatives: d_renin, d_aldosterone, d_ADH, d_ANP
    Also returns algebraic variables: angII, RSNA
    """
    if drug_effects is None:
        drug_effects = {}

    ACE_inhibition = drug_effects.get("ACE_inhibitor", 1.0)
    ARB_factor = drug_effects.get("ARB", 1.0)  # 1.0 = no block, 0.0 = full block

    # RSNA (algebraic — fast dynamics, treated as quasi-steady)
    RSNA = compute_RSNA(MAP, BV, params)

    # Ang II (quasi-steady due to very short half-life)
    angII = compute_angII(state["renin"], params, ACE_inhibition)

    # Effective Ang II (after ARB — receptor-level blockade)
    angII_effective = angII * ARB_factor

    # Renin secretion rate
    renin_secretion = compute_renin_secretion_rate(
        MAP, Na_to_DT, RSNA, angII, state["ANP"], params)

    # Derivatives
    d_renin = compute_renin_derivatives(state["renin"], renin_secretion, params)
    d_aldo = compute_aldosterone_derivatives(
        state["aldosterone"], angII_effective, state["ANP"], params)
    d_ADH = compute_ADH_derivatives(
        state["ADH"], plasma_osm, BV, angII_effective, params)
    d_ANP = compute_ANP_derivatives(state["ANP"], BV, params)

    return {
        "d_renin": d_renin,
        "d_aldosterone": d_aldo,
        "d_ADH": d_ADH,
        "d_ANP": d_ANP,
        "angII": angII,
        "angII_effective": angII_effective,
        "RSNA": RSNA,
        "renin_secretion": renin_secretion,
    }
