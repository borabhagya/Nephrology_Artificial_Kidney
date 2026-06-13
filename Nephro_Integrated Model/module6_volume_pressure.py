"""
module6_volume_pressure.py — Volume/Pressure Homeostasis Module

Computes the cardiovascular coupling that closes the whole-body feedback loop:

  Sodium balance → ECFV → Blood Volume → Cardiac Output → MAP → Kidney

Key equations:
  dNa_total/dt = Na_intake − Na_excreted
  dWater_total/dt = Water_intake − Urine_flow − Insensible_losses
  ECFV = f(Na_total, Water_total)
  BV = ECFV × BV_to_ECFV_ratio
  CO = f(BV) via Frank-Starling
  MAP = CO × SVR

This module represents the "long loop" that makes the kidney a blood pressure
regulator — the central insight of Guyton's model.

References:
  - Guyton et al. (1972) — original circulatory regulation model
  - Guyton & Hall, Ch. 19: Role of the Kidneys in Long-Term Control of
    Arterial Pressure
  - Karaaslan et al. (2005) — volume-pressure equations
"""

import numpy as np


def compute_plasma_sodium(total_Na, total_water):
    """
    Plasma sodium concentration (mEq/L).

    [Na] = total_body_Na / ECFV (in liters)

    This is a simplification: total body Na is largely extracellular,
    and plasma [Na] reflects the ratio of Na to water.
    """
    if total_water > 0:
        plasma_Na = total_Na / (total_water / 1000.0)  # mEq/L
    else:
        plasma_Na = 140.0
    return np.clip(plasma_Na, 100.0, 180.0)


def compute_plasma_osmolality(plasma_Na, params):
    """
    Plasma osmolality ≈ 2.1 × [Na] (simplified).

    Full formula: Osm = 2×[Na] + glucose/18 + BUN/2.8
    At normal glucose and BUN, ≈ 2.1 × [Na]
    """
    osm = params.volp.Na_to_osm_factor * plasma_Na
    return np.clip(osm, 220.0, 360.0)


def compute_blood_volume(ECFV, params):
    """
    Blood volume from extracellular fluid volume.

    BV = ECFV × BV_to_ECFV_ratio

    This ratio (~0.36) represents the fraction of ECF that is
    intravascular (plasma volume + RBC volume).
    """
    BV = ECFV * params.volp.BV_to_ECFV_ratio
    return max(BV, 2000.0)  # Physiological floor


def compute_cardiac_output(BV, angII, params):
    """
    Cardiac output from blood volume (simplified Frank-Starling).

    Increased BV → increased venous return → increased RAP → increased CO.
    Ang II increases systemic vascular resistance, indirectly modulating.

    CO = CO_basal + CO_slope × (RAP − RAP_threshold)
    RAP ≈ f(BV) via venous compliance

    Returns: CO (mL/min), RAP (mmHg)
    """
    vp = params.volp
    hemo = params.hemo

    # RAP from blood volume (simplified venous compliance model)
    # RAP = (BV − unstressed_volume) / venous_compliance
    # Calibrated so BV=5000 → RAP≈2 mmHg → CO≈5000 → MAP≈93
    unstressed_volume = 4800.0  # mL — volume at zero transmural pressure
    RAP = (BV - unstressed_volume) / hemo.venous_compliance
    RAP = np.clip(RAP, -2.0, 15.0)

    # Cardiac output from Frank-Starling
    CO = vp.CO_basal + vp.CO_slope * max(RAP - vp.RAP_threshold, 0.0)
    CO = np.clip(CO, 2000.0, 10000.0)

    return CO, RAP


def compute_MAP(CO, angII, RSNA, params):
    """
    Mean arterial pressure from cardiac output and systemic vascular resistance.

    MAP = CO × SVR

    SVR is modulated by:
      - Ang II (vasoconstriction)
      - RSNA (vasoconstriction)
    """
    hemo = params.hemo
    h = params.horm

    # SVR modulation
    angII_svr = 1.0 + h.angII_SVR_sensitivity * (angII - 1.0)
    rsna_svr = 1.0 + 0.03 * (RSNA - 1.0)

    SVR = hemo.SVR_normal * angII_svr * rsna_svr
    SVR = max(SVR, 0.005)

    MAP = CO * SVR  # CO in mL/min × SVR in mmHg/(mL/min) = mmHg
    MAP = np.clip(MAP, 40.0, 200.0)

    return MAP, SVR


def compute_volume_derivatives(Na_intake, water_intake, Na_excreted, urine_flow, params):
    """
    Compute time derivatives for sodium and water balance.

    dNa_total/dt = Na_intake − Na_excreted  (mEq/min)
    dWater_total/dt = Water_intake − Urine_flow − Insensible_losses  (mL/min)

    Parameters
    ----------
    Na_intake : float — Dietary sodium intake (mEq/min)
    water_intake : float — Total water intake (mL/min)
    Na_excreted : float — Urinary sodium excretion (mEq/min)
    urine_flow : float — Urine flow rate (mL/min)
    params : Parameters

    Returns
    -------
    dict with d_Na_total, d_water_total
    """
    insensible = params.tubW.insensible_loss / 1440.0  # Convert L/day → mL/min
    insensible *= 1000.0  # L to mL: 0.9 L/day = 0.625 mL/min

    d_Na_total = Na_intake - Na_excreted
    d_water_total = water_intake - urine_flow - insensible

    return {
        "d_Na_total": d_Na_total,
        "d_water_total": d_water_total,
    }


def compute_cardiovascular_state(total_Na, total_water, angII, RSNA, params):
    """
    Compute all cardiovascular variables from total body sodium and water.

    Returns
    -------
    dict with:
        plasma_Na, plasma_osm, ECFV, BV, CO, RAP, MAP, SVR
    """
    plasma_Na = compute_plasma_sodium(total_Na, total_water)
    plasma_osm = compute_plasma_osmolality(plasma_Na, params)

    ECFV = total_water  # ECFV ≈ total body water in this simplified model
    BV = compute_blood_volume(ECFV, params)

    CO, RAP = compute_cardiac_output(BV, angII, params)
    MAP, SVR = compute_MAP(CO, angII, RSNA, params)

    return {
        "plasma_Na": plasma_Na,
        "plasma_osm": plasma_osm,
        "ECFV": ECFV,
        "BV": BV,
        "CO": CO,
        "RAP": RAP,
        "MAP": MAP,
        "SVR": SVR,
    }
