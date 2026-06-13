"""
scenarios.py — Perturbation Scenarios for Model Validation

Defines standardized physiological and pharmacological perturbation scenarios
to test the integrated kidney model against known physiological responses.

Each scenario specifies:
  - Input functions (Na intake, water intake) as functions of time
  - Drug effect functions (if applicable)
  - Expected qualitative responses for validation
  - Simulation duration

References:
  - Guyton & Hall, Ch. 28–30 (renal responses to perturbations)
  - Karaaslan et al. (2005) — salt loading and RSNA scenarios
  - Hallow et al. (2014) — drug effect simulations
"""

import numpy as np
from parameters import Parameters


def _step_function(t, t_onset, baseline, perturbed):
    """Step from baseline to perturbed at t_onset."""
    return perturbed if t >= t_onset else baseline


# ── Scenario Definitions ──

def scenario_baseline(params=None):
    """
    Normal steady-state: constant normal Na and water intake.
    Verifies that the model converges to physiological equilibrium.

    Expected outcomes:
      - GFR ≈ 125 mL/min
      - MAP ≈ 93 mmHg
      - Urine flow ≈ 1 mL/min
      - Na excretion ≈ 150 mEq/day
      - All hormones at normalized 1.0
    """
    if params is None:
        params = Parameters()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    return {
        "name": "Baseline (Normal Steady State)",
        "params": params,
        "t_span": (0, 14400),  # 10 days
        "Na_intake_func": lambda t: Na_rate,
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": lambda t: {},
        "expected": {
            "GFR": (100, 150),
            "MAP": (85, 100),
            "urine_flow": (0.5, 2.0),
            "Na_excretion_daily": (100, 200),
            "plasma_Na": (135, 145),
        },
    }


def scenario_high_salt(params=None, salt_multiplier=3.0, t_onset=1440.0):
    """
    High salt intake: 3× normal Na intake starting at day 1.

    Expected responses (Guyton):
      - ↑ Na excretion to match new intake (pressure natriuresis)
      - ↓ Renin, ↓ Ang II, ↓ Aldosterone (RAAS suppression)
      - ↑ ANP
      - Modest ↑ ECFV, ↑ BV
      - Small ↑ MAP (salt-sensitive subjects may show more)
      - ↑ GFR (slightly, due to volume expansion)
      - New steady state within 2–5 days
    """
    if params is None:
        params = Parameters()

    Na_normal = params.tubNa.Na_intake_normal / 1440.0
    Na_high = Na_normal * salt_multiplier
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    return {
        "name": f"High Salt Intake ({salt_multiplier}× normal)",
        "params": params,
        "t_span": (0, 14400),
        "Na_intake_func": lambda t: _step_function(t, t_onset, Na_normal, Na_high),
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": lambda t: {},
        "expected": {
            "Na_excretion_daily_final": (350, 550),
            "MAP_increase": (0, 10),
            "renin_direction": "decrease",
            "aldosterone_direction": "decrease",
            "ANP_direction": "increase",
        },
    }


def scenario_low_salt(params=None, salt_multiplier=0.2, t_onset=1440.0):
    """
    Low salt intake: 0.2× normal (salt restriction).

    Expected responses:
      - ↓ Na excretion to match low intake
      - ↑ Renin, ↑ Ang II, ↑ Aldosterone (RAAS activation)
      - ↓ ANP
      - Modest ↓ ECFV, ↓ BV
      - Small ↓ MAP
    """
    if params is None:
        params = Parameters()

    Na_normal = params.tubNa.Na_intake_normal / 1440.0
    Na_low = Na_normal * salt_multiplier
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    return {
        "name": f"Low Salt Intake ({salt_multiplier}× normal)",
        "params": params,
        "t_span": (0, 14400),
        "Na_intake_func": lambda t: _step_function(t, t_onset, Na_normal, Na_low),
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": lambda t: {},
        "expected": {
            "renin_direction": "increase",
            "aldosterone_direction": "increase",
            "Na_excretion_daily_final": (20, 50),
        },
    }


def scenario_dehydration(params=None, water_multiplier=0.3, t_onset=1440.0):
    """
    Water deprivation: reduce water intake to 30% of normal.

    Expected responses (Guyton, StatPearls — Vasopressin):
      - ↑ ADH (↑ plasma osmolality triggers osmoreceptors)
      - ↓ Urine flow (concentrated urine)
      - ↑ Urine osmolality (up to 1200 mOsm/kg with max ADH)
      - Relatively preserved plasma Na (initially)
      - ↓ BV (moderate)
    """
    if params is None:
        params = Parameters()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_normal = params.tubW.water_intake_normal / 1440.0 * 1000.0
    W_low = W_normal * water_multiplier

    return {
        "name": f"Dehydration (water at {water_multiplier}× normal)",
        "params": params,
        "t_span": (0, 7200),  # 5 days
        "Na_intake_func": lambda t: Na_rate,
        "water_intake_func": lambda t: _step_function(t, t_onset, W_normal, W_low),
        "drug_effects_func": lambda t: {},
        "expected": {
            "ADH_direction": "increase",
            "urine_flow_direction": "decrease",
            "urine_osm_direction": "increase",
        },
    }


def scenario_volume_expansion(params=None, extra_water=500.0, t_onset=1440.0,
                               duration=60.0):
    """
    Acute isotonic saline infusion: +500 mL over 1 hour.

    Expected responses:
      - ↑ BV, ↑ ECFV (acute)
      - ↑ ANP (atrial stretch)
      - ↓ ADH, ↓ Renin
      - ↑ GFR (slightly)
      - ↑ Na excretion, ↑ urine flow
      - Return to baseline over hours
    """
    if params is None:
        params = Parameters()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_normal = params.tubW.water_intake_normal / 1440.0 * 1000.0
    bolus_rate = extra_water / duration  # mL/min during infusion

    def water_func(t):
        if t_onset <= t < t_onset + duration:
            return W_normal + bolus_rate
        return W_normal

    return {
        "name": "Acute Volume Expansion (500 mL saline)",
        "params": params,
        "t_span": (0, 4320),  # 3 days
        "Na_intake_func": lambda t: Na_rate + (
            (extra_water * 0.154 / 1000.0) / duration if t_onset <= t < t_onset + duration else 0
        ),  # isotonic saline has 154 mEq/L Na
        "water_intake_func": water_func,
        "drug_effects_func": lambda t: {},
        "expected": {
            "ANP_direction": "increase",
            "Na_excretion_direction": "increase",
        },
    }


def scenario_ACE_inhibitor(params=None, ACE_block=0.2, t_onset=1440.0):
    """
    ACE inhibitor: reduces ACE activity to 20% of normal.

    Expected responses (Hallow et al. 2014):
      - ↓ Ang II (primary effect)
      - ↓ Aldosterone (less Ang II stimulation)
      - ↑ Renin (loss of Ang II negative feedback)
      - ↓ MAP (vasodilation + natriuresis)
      - ↑ Na excretion (initially, then new steady state)
      - ↓ Efferent resistance → ↓ GFR (in some settings)
    """
    if params is None:
        params = Parameters()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    def drug_func(t):
        if t >= t_onset:
            return {"ACE_inhibitor": ACE_block}
        return {}

    return {
        "name": f"ACE Inhibitor (ACE at {ACE_block*100:.0f}%)",
        "params": params,
        "t_span": (0, 14400),
        "Na_intake_func": lambda t: Na_rate,
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": drug_func,
        "expected": {
            "angII_direction": "decrease",
            "aldosterone_direction": "decrease",
            "renin_direction": "increase",
            "MAP_direction": "decrease",
        },
    }


def scenario_loop_diuretic(params=None, block_fraction=0.5, t_onset=1440.0):
    """
    Loop diuretic (e.g., furosemide): blocks 50% of LoH Na reabsorption.

    Expected responses:
      - ↑ Na excretion (acute natriuresis)
      - ↑ Urine flow
      - ↓ BV, ↓ ECFV
      - ↑ Renin, ↑ Ang II, ↑ Aldosterone (compensatory RAAS activation)
      - ↓ MAP
    """
    if params is None:
        params = Parameters()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    def drug_func(t):
        if t >= t_onset:
            return {"loop_diuretic": 1.0 - block_fraction}
        return {}

    return {
        "name": f"Loop Diuretic ({block_fraction*100:.0f}% LoH block)",
        "params": params,
        "t_span": (0, 14400),
        "Na_intake_func": lambda t: Na_rate,
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": drug_func,
        "expected": {
            "Na_excretion_direction": "increase",
            "renin_direction": "increase",
            "MAP_direction": "decrease",
        },
    }


def scenario_renal_artery_stenosis(params=None, stenosis_factor=0.6, t_onset=1440.0):
    """
    Renal artery stenosis: reduces effective renal perfusion pressure.
    Simulated by increasing afferent resistance.

    Expected responses:
      - ↓ RBF, ↓ GFR (initially)
      - ↑ Renin, ↑ Ang II, ↑ Aldosterone (Goldblatt hypertension)
      - Na retention → ↑ ECFV, ↑ BV
      - ↑ MAP (renovascular hypertension)
    """
    if params is None:
        params = Parameters()

    # Simulate stenosis by increasing baseline afferent resistance
    params_modified = params.copy()
    # We modify the baseline resistance fraction to simulate stenosis

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    def Na_intake(t):
        return Na_rate

    def water_intake(t):
        return W_rate

    # Stenosis modeled as increased R_aa_fraction (fixed anatomical restriction)
    def drug_func(t):
        if t >= t_onset:
            return {"stenosis_R_aa_multiplier": 1.0 / stenosis_factor}
        return {}

    return {
        "name": f"Renal Artery Stenosis ({(1-stenosis_factor)*100:.0f}% reduction)",
        "params": params_modified,
        "t_span": (0, 14400),
        "Na_intake_func": Na_intake,
        "water_intake_func": water_intake,
        "drug_effects_func": drug_func,
        "expected": {
            "renin_direction": "increase",
            "angII_direction": "increase",
            "MAP_direction": "increase",
        },
    }


def scenario_sympathetic_activation(params=None, rsna_multiplier=2.0, t_onset=1440.0):
    """
    Sustained sympathetic activation: RSNA doubled.
    Simulated by shifting the RSNA baseline.

    Expected responses:
      - ↑ RVR, ↓ RBF
      - ↑ Renin (β1-mediated)
      - ↑ Proximal Na reabsorption
      - ↓ Na excretion → Na retention
      - ↑ BV, ↑ MAP
    """
    if params is None:
        params = Parameters()

    # Modify RSNA baseline
    params_modified = params.copy()

    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    def drug_func(t):
        if t >= t_onset:
            return {"rsna_offset": rsna_multiplier - 1.0}
        return {}

    return {
        "name": f"Sympathetic Activation (RSNA ×{rsna_multiplier})",
        "params": params_modified,
        "t_span": (0, 14400),
        "Na_intake_func": lambda t: Na_rate,
        "water_intake_func": lambda t: W_rate,
        "drug_effects_func": drug_func,
        "expected": {
            "renin_direction": "increase",
            "Na_excretion_direction": "decrease_then_normalize",
            "MAP_direction": "increase",
        },
    }


def get_all_scenarios(params=None):
    """Return a dict of all available scenarios."""
    return {
        "baseline": scenario_baseline(params),
        "high_salt": scenario_high_salt(params),
        "low_salt": scenario_low_salt(params),
        "dehydration": scenario_dehydration(params),
        "volume_expansion": scenario_volume_expansion(params),
        "ACE_inhibitor": scenario_ACE_inhibitor(params),
        "loop_diuretic": scenario_loop_diuretic(params),
        "renal_stenosis": scenario_renal_artery_stenosis(params),
        "sympathetic": scenario_sympathetic_activation(params),
    }
