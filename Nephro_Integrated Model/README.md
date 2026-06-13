# Integrated Kidney Model — A Feedback-Controlled Whole-Kidney Simulation

## Overview

A modular, ODE-based computational model of the kidney as an integrated physiological system, built on the Guyton/Karaaslan framework for circulatory and renal regulation. The model couples six functional modules through shared signals and feedback loops to simulate whole-body sodium balance, water balance, blood pressure regulation, and hormonal homeostasis.

## Module Architecture

```
┌───────────────────────────────────────────────────────────┐
│ Module 1: Renal Hemodynamics                              │
│   RBF, RPF, R_aa, R_ea, P_gc                             │
├───────────────────────────────────────────────────────────┤
│ Module 2: Glomerular Filtration                           │
│   GFR, FF, filtered loads (Starling equation)             │
├───────────────────────────────────────────────────────────┤
│ Module 3: Tubular Sodium Transport                        │
│   PT / LoH / DT / CD reabsorption, TGF signal            │
├───────────────────────────────────────────────────────────┤
│ Module 4: Tubular Water Handling                          │
│   Segment water reabsorption, ADH-regulated CD            │
├───────────────────────────────────────────────────────────┤
│ Module 5: Hormonal Regulation                             │
│   Renin, Ang II, Aldosterone, ADH, ANP, RSNA             │
├───────────────────────────────────────────────────────────┤
│ Module 6: Volume/Pressure Homeostasis                     │
│   Na balance, water balance, ECFV, BV, CO, MAP           │
└───────────────────────────────────────────────────────────┘
```

## State Variables (ODE)

| Variable       | Description                    | Units        |
|----------------|--------------------------------|--------------|
| total_Na       | Total body sodium              | mEq          |
| total_water    | Total body water (≈ ECFV)      | mL           |
| renin          | Plasma renin (normalized)      | dimensionless|
| aldosterone    | Plasma aldosterone (normalized)| dimensionless|
| ADH            | Plasma vasopressin (normalized)| dimensionless|
| ANP            | Plasma ANP (normalized)        | dimensionless|

## Running

```bash
# All scenarios
python run_simulation.py

# Single scenario
python run_simulation.py --scenario high_salt

# Sensitivity analysis
python sensitivity_analysis.py
```

## Available Scenarios

| Scenario              | Description                                    |
|-----------------------|------------------------------------------------|
| baseline              | Normal steady state verification               |
| high_salt             | 3× Na intake                                   |
| low_salt              | 0.2× Na intake (salt restriction)              |
| dehydration           | 30% normal water intake                        |
| volume_expansion      | 500 mL isotonic saline bolus                   |
| ACE_inhibitor         | ACE blocked to 20% activity                    |
| loop_diuretic         | 50% Loop of Henle Na reabsorption block        |
| renal_stenosis        | 40% renal perfusion reduction                  |
| sympathetic           | 2× baseline RSNA                               |

## Key References

1. Guyton AC et al. (1972) Circulatory regulation model
2. Karaaslan F et al. (2005) Ann Biomed Eng 33:1607–1630
3. Hallow KM et al. (2014) Am J Physiol Renal Physiol 307:F775–F784
4. Kutumova et al. (2021) Front Physiol 12:746300
5. Guyton & Hall, Textbook of Medical Physiology (14th ed.)
