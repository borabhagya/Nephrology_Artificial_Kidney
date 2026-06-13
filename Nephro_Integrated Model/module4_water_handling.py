"""
module4_water_handling.py — Tubular Water Handling Module

Computes water reabsorption along the nephron segments and final urine output.

Key principles:
  - Proximal tubule: iso-osmotic reabsorption (water follows sodium)
  - Descending limb of Henle: water reabsorbed into hypertonic medulla
  - Ascending limb: impermeable to water (diluting segment)
  - Distal tubule: relatively impermeable to water
  - Collecting duct: ADH-regulated water permeability
    - Without ADH: dilute urine (50–100 mOsm/kg)
    - With max ADH: concentrated urine (up to 1200 mOsm/kg)

References:
  - Guyton & Hall, Ch. 29: Urine Concentration and Dilution
  - Bankir et al. (2017) — Vasopressin physiology
"""

import numpy as np


def compute_PT_water_reabsorption(filtered_water, frac_Na_PT):
    """
    Proximal tubule water reabsorption.

    Iso-osmotic: water reabsorption fraction matches sodium reabsorption
    fraction due to high water permeability and solute coupling.

    Returns: (water_reabsorbed_PT, water_to_LoH)
    """
    water_reabs = filtered_water * frac_Na_PT
    water_to_LoH = filtered_water - water_reabs
    return water_reabs, water_to_LoH


def compute_LoH_water_reabsorption(water_delivered, params):
    """
    Loop of Henle water reabsorption.

    - Descending limb: permeable to water → reabsorption into hypertonic medulla
    - Ascending limb: impermeable to water (active NaCl reabsorption only)
    - Net: ~15% of water delivered to LoH is reabsorbed

    Returns: (water_reabsorbed_LoH, water_to_DT)
    """
    t = params.tubW
    water_reabs = water_delivered * t.frac_water_LoH_descending
    water_to_DT = water_delivered - water_reabs
    return water_reabs, water_to_DT


def compute_DT_water(water_delivered):
    """
    Distal tubule water handling.

    Early distal tubule (diluting segment) is relatively impermeable to water.
    Minimal water reabsorption in this segment.

    Returns: (water_reabsorbed_DT ≈ 0, water_to_CD)
    """
    return 0.0, water_delivered


def compute_CD_water_reabsorption(water_delivered, ADH, params):
    """
    Collecting duct water reabsorption — the key ADH-regulated segment.

    ADH increases aquaporin-2 insertion into apical membrane of principal cells,
    dramatically increasing water permeability.

    Fractional reabsorption modeled as sigmoidal function of ADH:
      frac_CD = frac_base + (frac_max − frac_base) × ADH_effect

    At ADH = 0 → dilute urine (minimal CD water reabsorption)
    At ADH = max → concentrated urine (maximal CD water reabsorption)

    Returns: (water_reabsorbed_CD, urine_flow_rate)
    """
    t = params.tubW

    # ADH effect: Michaelis-Menten saturation with steep response
    # Km = 0.05 so that at ADH=1.0 (normal), effect ≈ 0.95 → frac ≈ 0.97
    # At ADH=0 → effect=0 → frac=base (diabetes insipidus: dilute urine)
    # At ADH=2+ → effect≈1 → frac≈max (maximally concentrated urine)
    Km_ADH = 0.05
    ADH_effect = ADH / (ADH + Km_ADH)
    ADH_effect = np.clip(ADH_effect, 0.0, 1.0)

    frac_CD = t.frac_water_CD_base + (t.frac_water_CD_max - t.frac_water_CD_base) * ADH_effect
    frac_CD = np.clip(frac_CD, 0.0, 0.99)

    water_reabs = water_delivered * frac_CD
    urine_flow = water_delivered - water_reabs

    # Ensure minimum urine flow (obligatory solute excretion)
    urine_flow = max(urine_flow, 0.3)  # ~0.3 mL/min minimum (~430 mL/day)

    return water_reabs, urine_flow, frac_CD


def compute_urine_osmolality(Na_excreted, urine_flow, params):
    """
    Estimate urine osmolality from solute excretion rate and urine flow.

    Urine_osm ≈ (total solute excretion rate) / urine_flow
    Major solutes: Na (and accompanying anion), urea, K, etc.
    Simplified: total solute ≈ 2 × Na_excreted + urea_excretion

    Typical total solute excretion: ~600 mOsm/day = 0.417 mOsm/min
    """
    t = params.tubW

    # Estimate total osmolar excretion (mOsm/min)
    # Na contributes ~300 mOsm/day of the ~600 total at normal diet
    urea_excretion = 0.21  # mOsm/min (~300 mOsm/day, relatively constant)
    Na_osm = Na_excreted * 2.0  # Na + accompanying anion (Cl−)
    total_osm_excretion = Na_osm + urea_excretion

    if urine_flow > 0:
        urine_osm = total_osm_excretion / (urine_flow / 1000.0)  # mOsm/L
    else:
        urine_osm = t.urine_osm_max

    # Clamp to physiological range
    urine_osm = np.clip(urine_osm, t.urine_osm_min, t.urine_osm_max)

    return urine_osm


def compute_water_handling(filtered_water, frac_Na_PT, Na_excreted, ADH, params):
    """
    Full water handling computation.

    Parameters
    ----------
    filtered_water : float — Filtered water volume (mL/min) = GFR
    frac_Na_PT : float — Proximal tubule Na reabsorption fraction
    Na_excreted : float — Sodium excretion rate (mEq/min)
    ADH : float — Normalized ADH level
    params : Parameters

    Returns
    -------
    dict with water reabsorption and urine values
    """
    # Segment-by-segment
    water_reabs_PT, water_to_LoH = compute_PT_water_reabsorption(
        filtered_water, frac_Na_PT)

    water_reabs_LoH, water_to_DT = compute_LoH_water_reabsorption(
        water_to_LoH, params)

    water_reabs_DT, water_to_CD = compute_DT_water(water_to_DT)

    water_reabs_CD, urine_flow, frac_CD = compute_CD_water_reabsorption(
        water_to_CD, ADH, params)

    # Total water reabsorbed
    water_total_reabs = water_reabs_PT + water_reabs_LoH + water_reabs_DT + water_reabs_CD
    frac_total = water_total_reabs / filtered_water if filtered_water > 0 else 0.0

    # Urine osmolality
    urine_osm = compute_urine_osmolality(Na_excreted, urine_flow, params)

    # Free water clearance
    # CH2O = V − Cosm, where Cosm = Uosm × V / Posm
    posm = params.volp.plasma_osm_normal  # approximate, will be updated in module 6
    Cosm = urine_osm * (urine_flow / 1000.0) / posm * 1000.0 if posm > 0 else 0
    free_water_clearance = urine_flow - Cosm

    return {
        "water_reabs_PT": water_reabs_PT,
        "water_reabs_LoH": water_reabs_LoH,
        "water_reabs_DT": water_reabs_DT,
        "water_reabs_CD": water_reabs_CD,
        "water_to_LoH": water_to_LoH,
        "water_to_DT": water_to_DT,
        "water_to_CD": water_to_CD,
        "urine_flow": urine_flow,       # mL/min
        "urine_osm": urine_osm,         # mOsm/kg
        "water_total_reabs": water_total_reabs,
        "frac_total_reabs": frac_total,
        "frac_CD_water": frac_CD,
        "free_water_clearance": free_water_clearance,
    }
