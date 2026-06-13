"""
parameters.py — Research-calibrated physiological parameters for the integrated kidney model.

Parameter sources:
  - Guyton & Hall, Textbook of Medical Physiology (14th ed.)
  - Karaaslan et al., Ann Biomed Eng 33:1607–1630, 2005
  - Hallow et al., Am J Physiol Renal Physiol 307:F775–F784, 2014
  - Kutumova et al., Front Physiol 12:746300, 2021
  - StatPearls (NCBI Bookshelf) — GFR, RAAS, Vasopressin physiology
  - Bankir et al., J Intern Med 282:284–297, 2017 (vasopressin)

All values are for a reference 70 kg healthy adult male at steady state.
Units: mL, min, mmHg, mEq, mOsm, pg, ng unless otherwise noted.
"""

import dataclasses
from dataclasses import dataclass, field


@dataclass
class HemodynamicParams:
    """Renal and systemic hemodynamic parameters."""

    # --- Systemic cardiovascular ---
    MAP_normal: float = 93.0            # Mean arterial pressure (mmHg) — Guyton
    CO_normal: float = 5000.0           # Cardiac output (mL/min) — Guyton
    BV_normal: float = 5000.0           # Blood volume (mL)
    SVR_normal: float = 0.0186          # Systemic vascular resistance (mmHg·min/mL)
                                        # MAP / CO = 93/5000

    # --- Renal hemodynamics ---
    RBF_normal: float = 1100.0          # Renal blood flow (mL/min) — Guyton: ~22% of CO
    RPF_normal: float = 660.0           # Renal plasma flow (mL/min) — RBF × (1 − Hct)
    Hct: float = 0.40                   # Hematocrit
    RVR_normal: float = 0.0845          # Renal vascular resistance (mmHg·min/mL)
                                        # MAP / RBF ≈ 93/1100

    # --- Arteriolar resistances (fraction of total RVR) ---
    R_aa_fraction: float = 0.55         # Afferent arteriole fraction of RVR
    R_ea_fraction: float = 0.45         # Efferent arteriole fraction of RVR

    # --- Autoregulation ---
    autoregulation_lower: float = 80.0  # Lower limit of autoregulatory range (mmHg)
    autoregulation_upper: float = 170.0 # Upper limit (mmHg)
    myogenic_gain: float = 0.02         # Myogenic response gain (per mmHg)

    # --- Venous ---
    venous_compliance: float = 100.0    # Systemic venous compliance (mL/mmHg) — simplified
    RAP_normal: float = 2.0             # Right atrial pressure (mmHg)


@dataclass
class FiltrationParams:
    """Glomerular filtration parameters."""

    GFR_normal: float = 125.0           # Glomerular filtration rate (mL/min) — StatPearls
    Kf: float = 12.5                    # Ultrafiltration coefficient (mL/min/mmHg)
                                        # Guyton: ~12.5, ensures GFR = Kf × Pnet ≈ 125
    P_gc_normal: float = 60.0           # Glomerular capillary hydrostatic pressure (mmHg)
    P_bc: float = 18.0                  # Bowman's capsule pressure (mmHg) — Guyton
    pi_gc_normal: float = 32.0          # Average glomerular oncotic pressure (mmHg)
                                        # Entrance ~28, exit ~36, average ~32
    FF_normal: float = 0.19             # Filtration fraction — GFR/RPF ≈ 125/660

    # Oncotic pressure relationship (Landis-Pappenheimer)
    pi_coeff: float = 1.63              # π = pi_coeff × C_protein (mmHg per g/dL)
    plasma_protein: float = 7.0         # Plasma protein concentration (g/dL) — normal


@dataclass
class TubularSodiumParams:
    """Tubular sodium reabsorption parameters by segment."""

    # Fractional reabsorption at baseline (fraction of DELIVERED load to each segment)
    # Calibrated to match Guyton textbook: PT 67%, LoH 25%, DT 5%, CD 3% of FILTERED
    frac_reabs_PT: float = 0.67         # Proximal tubule — ~67% of filtered Na
    frac_reabs_LoH: float = 0.76        # Loop of Henle — 76% of delivered ≈ 25% of filtered
    frac_reabs_DT: float = 0.50         # Distal tubule — 50% of delivered ≈ 4% of filtered
    frac_reabs_CD: float = 0.85         # Collecting duct — 85% of delivered ≈ 3% of filtered
    # Net: ~99.4% reabsorbed → ~0.6% excreted ≈ 150 mEq/day at normal GFR

    # Filtered sodium load
    plasma_Na: float = 140.0            # Plasma sodium (mEq/L) — normal
    filtered_Na_normal: float = 25200.0 # Filtered Na (mEq/day) = GFR × P_Na × 1440/1000
                                        # = 125 × 140 × 1.44 ≈ 25,200 mEq/day

    # Hormonal modulation sensitivity coefficients
    # These represent the fractional change per unit normalized hormone change
    PT_angII_sensitivity: float = 0.10  # Ang II ↑ → ↑ PT reabsorption
    PT_rsna_sensitivity: float = 0.05   # RSNA ↑ → ↑ PT reabsorption
    PT_pressure_sensitivity: float = 0.003  # Pressure natriuresis in PT (per mmHg above normal)

    DT_aldosterone_sensitivity: float = 0.08  # Aldosterone ↑ → ↑ DT reabsorption
    CD_aldosterone_sensitivity: float = 0.12  # Aldosterone ↑ → ↑ CD reabsorption
    CD_anp_sensitivity: float = 0.10    # ANP ↑ → ↓ CD reabsorption

    # Normal sodium excretion (must match intake at steady state)
    Na_excretion_normal: float = 150.0  # mEq/day (~typical Western diet)
    Na_intake_normal: float = 150.0     # mEq/day


@dataclass
class TubularWaterParams:
    """Tubular water handling parameters."""

    # Fractional water reabsorption at baseline
    frac_water_PT: float = 0.67         # Iso-osmotic with Na in proximal tubule
    frac_water_LoH_descending: float = 0.15  # Descending limb — water reabsorbed
    frac_water_DT: float = 0.0          # Minimal in diluting segment
    frac_water_CD_base: float = 0.60    # Collecting duct baseline (without ADH — DI state)
    frac_water_CD_max: float = 0.993    # Max CD water reabsorption (full ADH, concentrated urine)

    # Urine parameters
    urine_flow_normal: float = 1.0      # mL/min (~1.5 L/day) — Guyton
    urine_osm_normal: float = 600.0     # mOsm/kg — moderately concentrated
    urine_osm_min: float = 50.0         # Maximum dilution (mOsm/kg)
    urine_osm_max: float = 1200.0       # Maximum concentration (mOsm/kg)

    # Water intake (calibrated to balance urine + insensible at steady state)
    water_intake_normal: float = 2.18    # L/day (mL/min ≈ 1.51)
    insensible_loss: float = 0.9         # L/day (skin + respiratory)

    # ADH sensitivity for CD water reabsorption
    ADH_water_sensitivity: float = 0.10  # fractional increase per pg/mL ADH


@dataclass
class HormonalParams:
    """Hormonal kinetics and stimulus-response parameters."""

    # --- Renin ---
    # Normal plasma renin activity: ~1.0 ng/mL/hr (supine, normal Na diet)
    # Half-life: 10–15 min in humans (Springer)
    renin_normal: float = 1.0           # Normalized to 1.0 (dimensionless)
    renin_half_life: float = 15.0       # min — normal human t1/2
    renin_secretion_rate_normal: float = 1.0  # Normalized
    renin_MAP_sensitivity: float = 0.011     # ↓ MAP → ↑ renin (inverse relationship)
    renin_macula_densa_sensitivity: float = 0.7  # ↓ MD Na → ↑ renin (inverse)
    renin_rsna_sensitivity: float = 0.5  # ↑ RSNA → ↑ renin (direct, β1)
    renin_angII_feedback: float = 0.3    # Ang II negative feedback on renin

    # --- Angiotensin II ---
    # Normal plasma: ~25 pg/mL; half-life <60 sec (StatPearls)
    angII_normal: float = 1.0           # Normalized
    angII_half_life: float = 0.5        # min (~30 sec, very short)
    angII_from_renin_gain: float = 1.0  # Proportional to renin (via ACE)

    # Ang II effects
    angII_efferent_sensitivity: float = 0.4   # ↑ Ang II → ↑ efferent resistance
    angII_afferent_sensitivity: float = 0.1   # Mild afferent constriction at high levels
    angII_SVR_sensitivity: float = 0.05       # Systemic vasoconstriction

    # --- Aldosterone ---
    # Normal plasma: ~5–15 ng/dL; half-life: ~20 min (hepatic clearance)
    aldo_normal: float = 1.0            # Normalized
    aldo_half_life: float = 20.0        # min
    aldo_from_angII_gain: float = 0.8   # Ang II stimulates aldosterone
    aldo_from_K_gain: float = 0.2       # Potassium effect (simplified, minor in this model)
    aldo_ANP_inhibition: float = 0.3    # ANP inhibits aldosterone secretion

    # --- ADH (Vasopressin / AVP) ---
    # Normal plasma: 1–3 pg/mL; half-life: ~3 min (Bankir 2017)
    # Osmotic threshold: 280 mOsm/kg (ScienceDirect/StatPearls)
    ADH_normal: float = 1.0             # Normalized
    ADH_half_life: float = 4.0          # min (conservative estimate, 3–5 min range)
    ADH_osm_threshold: float = 280.0    # mOsm/kg
    ADH_osm_sensitivity: float = 0.5    # Steep response above threshold
    ADH_volume_sensitivity: float = 0.2 # ↓ BV → ↑ ADH (requires >10-20% drop)
    ADH_angII_sensitivity: float = 0.1  # Ang II stimulates ADH release

    # --- ANP (Atrial Natriuretic Peptide) ---
    # Normal plasma: ~20–40 pg/mL; half-life: 2–5 min (StatPearls)
    ANP_normal: float = 1.0             # Normalized
    ANP_half_life: float = 3.0          # min
    ANP_volume_sensitivity: float = 0.8 # ↑ BV → ↑ ANP (atrial stretch)
    ANP_renin_inhibition: float = 0.3   # ANP inhibits renin secretion

    # --- RSNA (Renal Sympathetic Nerve Activity) ---
    RSNA_normal: float = 1.0            # Normalized
    RSNA_MAP_sensitivity: float = 0.012 # Baroreflex: ↓ MAP → ↑ RSNA
    RSNA_RAP_sensitivity: float = 0.03  # Cardiopulmonary reflex
    RSNA_time_constant: float = 10.0    # Response time (min)

    # Afferent arteriolar effects of RSNA
    RSNA_afferent_sensitivity: float = 0.2   # ↑ RSNA → ↑ R_aa
    RSNA_efferent_sensitivity: float = 0.1   # ↑ RSNA → ↑ R_ea


@dataclass
class TGFParams:
    """Tubuloglomerular feedback parameters."""

    # Sigmoidal TGF response
    TGF_midpoint: float = 1.0           # Normalized macula densa Na delivery at midpoint
    TGF_slope: float = 3.0              # Steepness of sigmoidal response
    TGF_max_constriction: float = 0.3   # Maximum fractional increase in R_aa from TGF
    TGF_max_dilation: float = 0.15      # Maximum fractional decrease in R_aa from TGF
    TGF_time_delay: float = 0.1         # Transit delay (min) — very fast


@dataclass
class VolumePressureParams:
    """Volume-pressure homeostasis / cardiovascular coupling."""

    # --- Extracellular fluid ---
    ECFV_normal: float = 14000.0        # Extracellular fluid volume (mL) — ~20% body weight
    total_body_Na_normal: float = 2000.0  # Total exchangeable Na (mEq) — approximate
    plasma_osm_normal: float = 290.0    # Plasma osmolality (mOsm/kg)

    # --- Volume distribution ---
    BV_to_ECFV_ratio: float = 0.357     # Blood volume / ECFV ≈ 5000/14000

    # --- Cardiac function (simplified Frank-Starling) ---
    # CO = CO_slope × (RAP − RAP_threshold) + CO_basal
    CO_slope: float = 500.0             # mL/min per mmHg RAP
    CO_basal: float = 4000.0            # Basal cardiac output (mL/min)
    RAP_threshold: float = 0.0          # mmHg

    # --- MAP from CO and SVR ---
    # MAP = CO × SVR  (this is the fundamental relation)

    # --- Sodium concentration and osmolality ---
    Na_to_osm_factor: float = 2.1       # Osm ≈ 2 × [Na] + glucose/18 + BUN/2.8
                                        # Simplified: Osm ≈ 2.1 × [Na] (at normal glucose/BUN)


@dataclass
class ModelConfig:
    """Simulation configuration."""

    dt: float = 0.1                     # Time step (min) — for ODE solver adaptive
    t_max: float = 14400.0              # Default simulation time (min) = 10 days
    solver_method: str = "Radau"        # Stiff-capable solver
    rtol: float = 1e-6                  # Relative tolerance
    atol: float = 1e-8                  # Absolute tolerance


@dataclass
class Parameters:
    """Master parameter container."""

    hemo: HemodynamicParams = field(default_factory=HemodynamicParams)
    filt: FiltrationParams = field(default_factory=FiltrationParams)
    tubNa: TubularSodiumParams = field(default_factory=TubularSodiumParams)
    tubW: TubularWaterParams = field(default_factory=TubularWaterParams)
    horm: HormonalParams = field(default_factory=HormonalParams)
    tgf: TGFParams = field(default_factory=TGFParams)
    volp: VolumePressureParams = field(default_factory=VolumePressureParams)
    config: ModelConfig = field(default_factory=ModelConfig)

    def copy(self):
        """Return a deep copy of all parameters."""
        import copy
        return copy.deepcopy(self)
