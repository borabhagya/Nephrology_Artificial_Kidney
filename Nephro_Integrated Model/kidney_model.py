"""
kidney_model.py — Integrated Whole-Kidney Model

Connects all six modules into a single coupled ODE system:
  Module 1: Renal Hemodynamics
  Module 2: Glomerular Filtration
  Module 3: Tubular Sodium Transport
  Module 4: Tubular Water Handling
  Module 5: Hormonal Regulation (RAAS, ADH, ANP, RSNA)
  Module 6: Volume/Pressure Homeostasis

State vector y = [total_Na, total_water, renin, aldosterone, ADH, ANP]
  - total_Na: total body sodium (mEq)
  - total_water: total body water / ECFV (mL)
  - renin: normalized plasma renin (dimensionless)
  - aldosterone: normalized plasma aldosterone (dimensionless)
  - ADH: normalized plasma ADH (dimensionless)
  - ANP: normalized plasma ANP (dimensionless)

Algebraic (computed at each time step):
  - RSNA, angII (quasi-steady-state)
  - MAP, BV, CO, plasma_Na, plasma_osm
  - GFR, RBF, R_aa, R_ea
  - Segmental Na reabsorption, Na excretion
  - Urine flow, urine osmolality
  - TGF signal

References:
  - Guyton et al. (1972) — Systems model of circulatory regulation
  - Karaaslan et al. (2005) — Long-term RSNA-renal model
  - Hallow et al. (2014) — Extended Guyton/Karaaslan model
"""

import numpy as np
from scipy.integrate import solve_ivp
from parameters import Parameters
from module1_hemodynamics import compute_renal_hemodynamics
from module2_filtration import compute_GFR
from module3_sodium_transport import compute_sodium_transport
from module4_water_handling import compute_water_handling
from module5_hormonal import (compute_RSNA, compute_angII,
                               compute_hormonal_derivatives)
from module6_volume_pressure import (compute_cardiovascular_state,
                                      compute_volume_derivatives)


# ── State vector indexing ──
STATE_NAMES = ["total_Na", "total_water", "renin", "aldosterone", "ADH", "ANP"]
IDX = {name: i for i, name in enumerate(STATE_NAMES)}
N_STATES = len(STATE_NAMES)


def get_initial_state(params):
    """
    Return the initial state vector for a healthy adult at steady state.
    """
    vp = params.volp
    h = params.horm

    y0 = np.zeros(N_STATES)
    y0[IDX["total_Na"]] = vp.total_body_Na_normal       # mEq
    y0[IDX["total_water"]] = vp.ECFV_normal              # mL
    y0[IDX["renin"]] = h.renin_normal                    # normalized
    y0[IDX["aldosterone"]] = h.aldo_normal               # normalized
    y0[IDX["ADH"]] = h.ADH_normal                        # normalized
    y0[IDX["ANP"]] = h.ANP_normal                        # normalized

    return y0


def compute_all_algebraic(y, params, drug_effects=None):
    """
    Given the current state vector, compute all algebraic (non-state) variables.

    This function evaluates all modules in the correct dependency order and
    returns a comprehensive snapshot of the system.

    Parameters
    ----------
    y : array-like — state vector
    params : Parameters
    drug_effects : dict or None — drug modulation factors

    Returns
    -------
    dict — all computed variables
    """
    if drug_effects is None:
        drug_effects = {}

    total_Na = y[IDX["total_Na"]]
    total_water = y[IDX["total_water"]]
    renin = max(y[IDX["renin"]], 0.01)
    aldosterone = max(y[IDX["aldosterone"]], 0.01)
    ADH = max(y[IDX["ADH"]], 0.01)
    ANP = max(y[IDX["ANP"]], 0.01)

    # ── Step 1: Cardiovascular state (Module 6 — partial) ──
    # First pass: compute RSNA and AngII with current state
    ACE_inhibition = drug_effects.get("ACE_inhibitor", 1.0)
    ARB_factor = drug_effects.get("ARB", 1.0)

    # Preliminary MAP estimate for RSNA (use previous or normal)
    # We'll iterate to self-consistency
    angII = compute_angII(renin, params, ACE_inhibition) * ARB_factor
    RSNA_prelim = compute_RSNA(params.hemo.MAP_normal, params.hemo.BV_normal, params)

    cv = compute_cardiovascular_state(total_Na, total_water, angII, RSNA_prelim, params)
    MAP = cv["MAP"]
    BV = cv["BV"]

    # Refine RSNA with actual MAP and BV
    RSNA = compute_RSNA(MAP, BV, params)

    # Recompute MAP with refined RSNA
    cv = compute_cardiovascular_state(total_Na, total_water, angII, RSNA, params)
    MAP = cv["MAP"]
    BV = cv["BV"]

    # ── Step 2: Initial TGF (use previous value or zero) ──
    TGF_signal = 0.0

    # ── Step 3: Renal hemodynamics (Module 1) ──
    hemo = compute_renal_hemodynamics(MAP, RSNA, angII, TGF_signal, params)

    # ── Step 4: Glomerular filtration (Module 2) ──
    filt = compute_GFR(hemo["P_gc"], hemo["RPF"], params)

    # ── Step 5: Sodium transport (Module 3) ──
    sodium = compute_sodium_transport(
        filt["filtered_Na"], MAP, angII, RSNA, aldosterone, ANP, params,
        drug_effects)

    # ── Step 6: Update TGF and re-iterate hemodynamics ──
    TGF_signal = sodium["TGF_signal"]
    hemo = compute_renal_hemodynamics(MAP, RSNA, angII, TGF_signal, params)
    filt = compute_GFR(hemo["P_gc"], hemo["RPF"], params)
    sodium = compute_sodium_transport(
        filt["filtered_Na"], MAP, angII, RSNA, aldosterone, ANP, params,
        drug_effects)

    # ── Step 7: Water handling (Module 4) ──
    water = compute_water_handling(
        filt["filtered_water"], sodium["frac_PT"],
        sodium["Na_excreted"], ADH, params)

    # ── Step 8: Plasma osmolality (for ADH feedback) ──
    plasma_Na = cv["plasma_Na"]
    plasma_osm = cv["plasma_osm"]

    # ── Compile all results ──
    results = {
        # State variables
        "total_Na": total_Na,
        "total_water": total_water,
        "renin": renin,
        "aldosterone": aldosterone,
        "ADH": ADH,
        "ANP": ANP,

        # Algebraic: hormonal
        "angII": angII,
        "RSNA": RSNA,

        # Cardiovascular
        "MAP": MAP,
        "BV": BV,
        "CO": cv["CO"],
        "RAP": cv["RAP"],
        "SVR": cv["SVR"],
        "ECFV": cv["ECFV"],
        "plasma_Na": plasma_Na,
        "plasma_osm": plasma_osm,

        # Hemodynamic
        "RBF": hemo["RBF"],
        "RPF": hemo["RPF"],
        "RVR": hemo["RVR"],
        "R_aa": hemo["R_aa"],
        "R_ea": hemo["R_ea"],
        "P_gc": hemo["P_gc"],

        # Filtration
        "GFR": filt["GFR"],
        "FF": filt["FF"],
        "P_net": filt["P_net"],
        "filtered_Na": filt["filtered_Na"],

        # Sodium transport
        "Na_excreted": sodium["Na_excreted"],
        "Na_excretion_daily": sodium["Na_excreted"] * 1440.0,  # mEq/day
        "frac_total_reabs": sodium["frac_total_reabs"],
        "TGF_signal": sodium["TGF_signal"],
        "Na_to_DT": sodium["Na_to_DT"],

        # Water handling
        "urine_flow": water["urine_flow"],
        "urine_flow_daily": water["urine_flow"] * 1.44,  # L/day
        "urine_osm": water["urine_osm"],
        "free_water_clearance": water["free_water_clearance"],
    }

    return results


def system_derivatives(t, y, params, Na_intake_func, water_intake_func,
                       drug_effects_func=None):
    """
    ODE right-hand side for the coupled kidney model.

    dy/dt = f(t, y, inputs)

    Parameters
    ----------
    t : float — current time (min)
    y : array — state vector
    params : Parameters
    Na_intake_func : callable(t) → float — Na intake rate (mEq/min)
    water_intake_func : callable(t) → float — water intake rate (mL/min)
    drug_effects_func : callable(t) → dict or None

    Returns
    -------
    dydt : array — time derivatives of state vector
    """
    # Ensure non-negative states
    y = np.maximum(y, 0.01)

    # Get inputs at current time
    Na_intake = Na_intake_func(t)
    water_intake = water_intake_func(t)
    drug_effects = drug_effects_func(t) if drug_effects_func else {}

    # Compute all algebraic variables
    alg = compute_all_algebraic(y, params, drug_effects)

    # Hormone state dict for Module 5
    hormone_state = {
        "renin": y[IDX["renin"]],
        "aldosterone": y[IDX["aldosterone"]],
        "ADH": y[IDX["ADH"]],
        "ANP": y[IDX["ANP"]],
    }

    # Hormonal derivatives (Module 5)
    horm_derivs = compute_hormonal_derivatives(
        hormone_state, alg["MAP"], alg["BV"], alg["plasma_osm"],
        alg["Na_to_DT"], params, drug_effects)

    # Volume derivatives (Module 6)
    vol_derivs = compute_volume_derivatives(
        Na_intake, water_intake, alg["Na_excreted"], alg["urine_flow"], params)

    # Assemble derivative vector
    dydt = np.zeros(N_STATES)
    dydt[IDX["total_Na"]] = vol_derivs["d_Na_total"]
    dydt[IDX["total_water"]] = vol_derivs["d_water_total"]
    dydt[IDX["renin"]] = horm_derivs["d_renin"]
    dydt[IDX["aldosterone"]] = horm_derivs["d_aldosterone"]
    dydt[IDX["ADH"]] = horm_derivs["d_ADH"]
    dydt[IDX["ANP"]] = horm_derivs["d_ANP"]

    return dydt


def run_simulation(params=None, t_span=None, y0=None,
                   Na_intake_func=None, water_intake_func=None,
                   drug_effects_func=None,
                   t_eval=None, max_step=10.0):
    """
    Run the integrated kidney model simulation.

    Parameters
    ----------
    params : Parameters (default: normal healthy adult)
    t_span : tuple (t_start, t_end) in minutes
    y0 : array — initial state (default: healthy steady state)
    Na_intake_func : callable(t) → float — Na intake (mEq/min)
    water_intake_func : callable(t) → float — water intake (mL/min)
    drug_effects_func : callable(t) → dict
    t_eval : array — time points for output
    max_step : float — max ODE solver step size (min)

    Returns
    -------
    solution : ODE solution object
    """
    if params is None:
        params = Parameters()

    if y0 is None:
        y0 = get_initial_state(params)

    if t_span is None:
        t_span = (0, params.config.t_max)

    # Default intake functions (steady state)
    if Na_intake_func is None:
        Na_rate = params.tubNa.Na_intake_normal / 1440.0  # mEq/day → mEq/min
        Na_intake_func = lambda t: Na_rate

    if water_intake_func is None:
        W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0  # L/day → mL/min
        water_intake_func = lambda t: W_rate

    if drug_effects_func is None:
        drug_effects_func = lambda t: {}

    # Solve ODE system
    sol = solve_ivp(
        fun=lambda t, y: system_derivatives(
            t, y, params, Na_intake_func, water_intake_func, drug_effects_func),
        t_span=t_span,
        y0=y0,
        method=params.config.solver_method,
        rtol=params.config.rtol,
        atol=params.config.atol,
        t_eval=t_eval,
        max_step=max_step,
        dense_output=True,
    )

    return sol


def extract_timeseries(sol, params, drug_effects_func=None):
    """
    Post-process ODE solution to extract all variables at each time point.

    Returns
    -------
    dict of arrays — each key is a variable name, value is array over time
    """
    if drug_effects_func is None:
        drug_effects_func = lambda t: {}

    n_times = len(sol.t)
    results = {}

    for i in range(n_times):
        y = sol.y[:, i]
        t = sol.t[i]
        drug_effects = drug_effects_func(t)
        alg = compute_all_algebraic(y, params, drug_effects)

        if i == 0:
            for key in alg:
                results[key] = np.zeros(n_times)
            results["time"] = np.zeros(n_times)
            results["time_days"] = np.zeros(n_times)

        results["time"][i] = t
        results["time_days"][i] = t / 1440.0
        for key, val in alg.items():
            results[key][i] = val

    return results
