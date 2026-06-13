"""
sensitivity_analysis.py — Parameter Sensitivity Analysis

Performs local (one-at-a-time) sensitivity analysis to identify which
parameters most strongly influence key model outputs.

Methodology:
  - Perturb each parameter by ±10% from baseline
  - Run to steady state
  - Compute normalized sensitivity coefficient:
    S_ij = (ΔOutput_j / Output_j) / (ΔParam_i / Param_i)

References:
  - Saltelli et al. (2008) — Global Sensitivity Analysis
  - Hester et al. (2011) — HumMod sensitivity analysis
  - PLOS Comp Biol (2012) — Virtual patients and sensitivity of Guyton model
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

from parameters import Parameters
from kidney_model import (run_simulation, extract_timeseries,
                          get_initial_state, compute_all_algebraic)


# Parameters to perturb and their access paths
PARAM_SPECS = [
    ("Kf (Ultrafilt Coeff)", "filt", "Kf"),
    ("R_aa fraction", "hemo", "R_aa_fraction"),
    ("Frac reabs PT", "tubNa", "frac_reabs_PT"),
    ("Frac reabs LoH", "tubNa", "frac_reabs_LoH"),
    ("CD aldo sensitivity", "tubNa", "CD_aldosterone_sensitivity"),
    ("CD ANP sensitivity", "tubNa", "CD_anp_sensitivity"),
    ("PT AngII sensitivity", "tubNa", "PT_angII_sensitivity"),
    ("Pressure natriuresis", "tubNa", "PT_pressure_sensitivity"),
    ("Renin t1/2", "horm", "renin_half_life"),
    ("AngII efferent sens", "horm", "angII_efferent_sensitivity"),
    ("Aldo from AngII gain", "horm", "aldo_from_angII_gain"),
    ("ADH osm sensitivity", "horm", "ADH_osm_sensitivity"),
    ("ANP volume sensitivity", "horm", "ANP_volume_sensitivity"),
    ("RSNA MAP sensitivity", "horm", "RSNA_MAP_sensitivity"),
    ("BV/ECFV ratio", "volp", "BV_to_ECFV_ratio"),
    ("CO slope", "volp", "CO_slope"),
    ("TGF slope", "tgf", "TGF_slope"),
    ("ADH water sens", "tubW", "ADH_water_sensitivity"),
    ("Myogenic gain", "hemo", "myogenic_gain"),
    ("SVR normal", "hemo", "SVR_normal"),
]

# Output variables to measure
OUTPUT_VARS = ["MAP", "GFR", "Na_excretion_daily", "urine_flow",
               "renin", "aldosterone", "ADH", "BV"]


def get_param_value(params, module, attr):
    """Get a parameter value by module name and attribute."""
    return getattr(getattr(params, module), attr)


def set_param_value(params, module, attr, value):
    """Set a parameter value by module name and attribute."""
    setattr(getattr(params, module), attr, value)


def run_steady_state(params, t_days=10):
    """Run model to steady state and return final outputs."""
    t_span = (0, t_days * 1440)
    Na_rate = params.tubNa.Na_intake_normal / 1440.0
    W_rate = params.tubW.water_intake_normal / 1440.0 * 1000.0

    t_eval = np.array([t_span[1] - 100, t_span[1]])

    sol = run_simulation(
        params=params,
        t_span=t_span,
        Na_intake_func=lambda t: Na_rate,
        water_intake_func=lambda t: W_rate,
        t_eval=t_eval,
        max_step=10.0,
    )

    ts = extract_timeseries(sol, params)
    return {var: ts[var][-1] for var in OUTPUT_VARS}


def compute_sensitivity_matrix(perturbation=0.10):
    """
    Compute normalized sensitivity matrix.

    S[i, j] = (ΔY_j / Y_j_baseline) / (ΔP_i / P_i_baseline)

    Parameters
    ----------
    perturbation : float — fractional perturbation (default 10%)

    Returns
    -------
    S : array (n_params × n_outputs)
    param_names : list of str
    output_names : list of str
    baseline_values : dict
    """
    print("Computing sensitivity matrix...")
    print(f"  Perturbation: ±{perturbation*100:.0f}%")
    print(f"  Parameters: {len(PARAM_SPECS)}")
    print(f"  Outputs: {len(OUTPUT_VARS)}")

    # Baseline run
    params_base = Parameters()
    print("\n  Running baseline...")
    baseline = run_steady_state(params_base)
    print(f"  Baseline MAP={baseline['MAP']:.1f}, GFR={baseline['GFR']:.1f}")

    n_params = len(PARAM_SPECS)
    n_outputs = len(OUTPUT_VARS)
    S = np.zeros((n_params, n_outputs))

    for i, (pname, module, attr) in enumerate(PARAM_SPECS):
        print(f"  [{i+1}/{n_params}] Perturbing: {pname}...", end="", flush=True)

        params_up = params_base.copy()
        params_dn = params_base.copy()

        base_val = get_param_value(params_base, module, attr)
        val_up = base_val * (1.0 + perturbation)
        val_dn = base_val * (1.0 - perturbation)

        set_param_value(params_up, module, attr, val_up)
        set_param_value(params_dn, module, attr, val_dn)

        try:
            out_up = run_steady_state(params_up, t_days=8)
            out_dn = run_steady_state(params_dn, t_days=8)

            for j, var in enumerate(OUTPUT_VARS):
                y_base = baseline[var]
                y_up = out_up[var]
                y_dn = out_dn[var]

                if abs(y_base) > 1e-10 and abs(base_val) > 1e-10:
                    # Central difference sensitivity
                    dy = (y_up - y_dn) / (2.0 * perturbation * y_base)
                    S[i, j] = dy
                else:
                    S[i, j] = 0.0

            print(f" done (MAP sens={S[i, 0]:.3f})")
        except Exception as e:
            print(f" FAILED: {e}")
            S[i, :] = 0.0

    return S, [p[0] for p in PARAM_SPECS], OUTPUT_VARS, baseline


def plot_sensitivity_heatmap(S, param_names, output_names, output_dir):
    """Plot sensitivity matrix as a heatmap."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Clip for visualization
    S_clip = np.clip(S, -5, 5)

    im = ax.imshow(S_clip, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(output_names)))
    ax.set_xticklabels(output_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(param_names)))
    ax.set_yticklabels(param_names, fontsize=9)

    # Add text annotations
    for i in range(len(param_names)):
        for j in range(len(output_names)):
            val = S[i, j]
            color = "white" if abs(val) > 1.0 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=7, color=color)

    ax.set_title("Parameter Sensitivity Matrix\n(Normalized: ΔOutput/Output per ΔParam/Param)",
                 fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax, label="Sensitivity Coefficient", shrink=0.8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sensitivity_heatmap.png"))
    plt.close()
    print(f"  Heatmap saved: sensitivity_heatmap.png")


def plot_tornado(S, param_names, output_names, output_dir):
    """Plot tornado diagrams for key outputs."""
    key_outputs = ["MAP", "GFR", "Na_excretion_daily"]

    fig, axes = plt.subplots(1, len(key_outputs), figsize=(15, 8))

    for k, var in enumerate(key_outputs):
        ax = axes[k]
        j = output_names.index(var)
        sensitivities = S[:, j]

        # Sort by absolute sensitivity
        order = np.argsort(np.abs(sensitivities))[::-1][:15]  # Top 15

        names = [param_names[i] for i in order]
        values = [sensitivities[i] for i in order]

        colors = ["#DC2626" if v > 0 else "#2563EB" for v in values]

        y_pos = np.arange(len(names))
        ax.barh(y_pos, values, color=colors, alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("Sensitivity Coefficient")
        ax.set_title(f"Tornado: {var}", fontweight="bold")
        ax.axvline(0, color="black", linewidth=0.5)
        ax.invert_yaxis()

    plt.suptitle("Parameter Sensitivity — Tornado Diagrams",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sensitivity_tornado.png"))
    plt.close()
    print(f"  Tornado plots saved: sensitivity_tornado.png")


def main():
    output_dir = "/home/claude/kidney_model/output"
    os.makedirs(output_dir, exist_ok=True)

    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   PARAMETER SENSITIVITY ANALYSIS                           ║")
    print("║   Local One-At-A-Time · Normalized Coefficients            ║")
    print("╚══════════════════════════════════════════════════════════════╝")

    S, param_names, output_names, baseline = compute_sensitivity_matrix(perturbation=0.10)

    print("\n  Generating plots...")
    plot_sensitivity_heatmap(S, param_names, output_names, output_dir)
    plot_tornado(S, param_names, output_names, output_dir)

    # Print top sensitivities for MAP
    print("\n  TOP 5 Parameters for MAP Sensitivity:")
    j_map = output_names.index("MAP")
    order = np.argsort(np.abs(S[:, j_map]))[::-1][:5]
    for rank, i in enumerate(order):
        print(f"    {rank+1}. {param_names[i]:30s}  S = {S[i, j_map]:+.3f}")

    print("\n  TOP 5 Parameters for GFR Sensitivity:")
    j_gfr = output_names.index("GFR")
    order = np.argsort(np.abs(S[:, j_gfr]))[::-1][:5]
    for rank, i in enumerate(order):
        print(f"    {rank+1}. {param_names[i]:30s}  S = {S[i, j_gfr]:+.3f}")

    print(f"\n✓ Sensitivity analysis complete. Outputs in: {output_dir}")


if __name__ == "__main__":
    main()
