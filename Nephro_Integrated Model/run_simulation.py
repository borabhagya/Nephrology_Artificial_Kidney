"""
run_simulation.py — Main entry point for running and visualizing kidney model simulations.

Usage:
  python run_simulation.py                    # Run all scenarios
  python run_simulation.py --scenario baseline  # Run specific scenario
  python run_simulation.py --list             # List available scenarios

Generates publication-quality validation plots for each scenario.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from parameters import Parameters
from kidney_model import (run_simulation, extract_timeseries,
                          get_initial_state, compute_all_algebraic)
from scenarios import get_all_scenarios


# ── Plot styling ──
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
})

COLORS = {
    "primary": "#2563EB",
    "secondary": "#DC2626",
    "tertiary": "#16A34A",
    "quaternary": "#D97706",
    "quinary": "#7C3AED",
    "senary": "#DB2777",
}


def print_steady_state_report(ts, label=""):
    """Print key steady-state values from the last time point."""
    idx = -1  # last time point
    print(f"\n{'='*60}")
    print(f"  STEADY-STATE REPORT: {label}")
    print(f"{'='*60}")
    print(f"  Time:              {ts['time_days'][idx]:.1f} days")
    print(f"")
    print(f"  --- Cardiovascular ---")
    print(f"  MAP:               {ts['MAP'][idx]:.1f} mmHg")
    print(f"  CO:                {ts['CO'][idx]:.0f} mL/min")
    print(f"  Blood Volume:      {ts['BV'][idx]:.0f} mL")
    print(f"  ECFV:              {ts['ECFV'][idx]:.0f} mL")
    print(f"  Plasma [Na]:       {ts['plasma_Na'][idx]:.1f} mEq/L")
    print(f"  Plasma Osm:        {ts['plasma_osm'][idx]:.0f} mOsm/kg")
    print(f"")
    print(f"  --- Renal Hemodynamics ---")
    print(f"  RBF:               {ts['RBF'][idx]:.0f} mL/min")
    print(f"  GFR:               {ts['GFR'][idx]:.1f} mL/min")
    print(f"  FF:                {ts['FF'][idx]:.3f}")
    print(f"")
    print(f"  --- Excretion ---")
    print(f"  Na Excretion:      {ts['Na_excretion_daily'][idx]:.1f} mEq/day")
    print(f"  Urine Flow:        {ts['urine_flow'][idx]:.2f} mL/min ({ts['urine_flow_daily'][idx]:.2f} L/day)")
    print(f"  Urine Osm:         {ts['urine_osm'][idx]:.0f} mOsm/kg")
    print(f"")
    print(f"  --- Hormones (normalized) ---")
    print(f"  Renin:             {ts['renin'][idx]:.3f}")
    print(f"  Ang II:            {ts['angII'][idx]:.3f}")
    print(f"  Aldosterone:       {ts['aldosterone'][idx]:.3f}")
    print(f"  ADH:               {ts['ADH'][idx]:.3f}")
    print(f"  ANP:               {ts['ANP'][idx]:.3f}")
    print(f"  RSNA:              {ts['RSNA'][idx]:.3f}")
    print(f"{'='*60}")


def plot_scenario(ts, scenario_info, output_dir):
    """Generate comprehensive validation plots for a scenario."""
    name = scenario_info["name"]
    safe_name = name.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
    safe_name = safe_name.replace("×", "x").replace("%", "pct")

    t_days = ts["time_days"]

    fig = plt.figure(figsize=(16, 20))
    fig.suptitle(f"Integrated Kidney Model — {name}", fontsize=14, fontweight="bold", y=0.98)

    gs = GridSpec(5, 3, figure=fig, hspace=0.4, wspace=0.35)

    # Row 1: Cardiovascular
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t_days, ts["MAP"], color=COLORS["primary"], linewidth=1.5)
    ax1.set_ylabel("MAP (mmHg)")
    ax1.set_title("Mean Arterial Pressure")
    ax1.axhline(93, color="gray", linestyle="--", alpha=0.5, label="Normal")
    ax1.legend(fontsize=8)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(t_days, ts["BV"], color=COLORS["secondary"], linewidth=1.5)
    ax2.set_ylabel("BV (mL)")
    ax2.set_title("Blood Volume")
    ax2.axhline(5000, color="gray", linestyle="--", alpha=0.5)

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(t_days, ts["plasma_Na"], color=COLORS["tertiary"], linewidth=1.5)
    ax3.set_ylabel("[Na] (mEq/L)")
    ax3.set_title("Plasma Sodium")
    ax3.axhline(140, color="gray", linestyle="--", alpha=0.5)

    # Row 2: Renal hemodynamics and filtration
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(t_days, ts["RBF"], color=COLORS["primary"], linewidth=1.5, label="RBF")
    ax4.set_ylabel("Flow (mL/min)")
    ax4.set_title("Renal Blood Flow")
    ax4.axhline(1100, color="gray", linestyle="--", alpha=0.5)

    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(t_days, ts["GFR"], color=COLORS["secondary"], linewidth=1.5)
    ax5.set_ylabel("GFR (mL/min)")
    ax5.set_title("Glomerular Filtration Rate")
    ax5.axhline(125, color="gray", linestyle="--", alpha=0.5)

    ax6 = fig.add_subplot(gs[1, 2])
    ax6.plot(t_days, ts["FF"], color=COLORS["tertiary"], linewidth=1.5)
    ax6.set_ylabel("FF")
    ax6.set_title("Filtration Fraction")
    ax6.axhline(0.19, color="gray", linestyle="--", alpha=0.5)

    # Row 3: Sodium and water excretion
    ax7 = fig.add_subplot(gs[2, 0])
    ax7.plot(t_days, ts["Na_excretion_daily"], color=COLORS["quaternary"], linewidth=1.5)
    ax7.set_ylabel("Na Excretion (mEq/day)")
    ax7.set_title("Sodium Excretion")
    ax7.axhline(150, color="gray", linestyle="--", alpha=0.5)

    ax8 = fig.add_subplot(gs[2, 1])
    ax8.plot(t_days, ts["urine_flow"], color=COLORS["quinary"], linewidth=1.5)
    ax8.set_ylabel("Urine Flow (mL/min)")
    ax8.set_title("Urine Flow Rate")
    ax8.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    ax9 = fig.add_subplot(gs[2, 2])
    ax9.plot(t_days, ts["urine_osm"], color=COLORS["senary"], linewidth=1.5)
    ax9.set_ylabel("Urine Osm (mOsm/kg)")
    ax9.set_title("Urine Osmolality")
    ax9.axhline(600, color="gray", linestyle="--", alpha=0.5)

    # Row 4: Hormones (RAAS)
    ax10 = fig.add_subplot(gs[3, 0])
    ax10.plot(t_days, ts["renin"], color=COLORS["primary"], linewidth=1.5, label="Renin")
    ax10.plot(t_days, ts["angII"], color=COLORS["secondary"], linewidth=1.5, label="Ang II")
    ax10.set_ylabel("Normalized Level")
    ax10.set_title("RAAS — Renin & Ang II")
    ax10.legend(fontsize=8)
    ax10.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    ax11 = fig.add_subplot(gs[3, 1])
    ax11.plot(t_days, ts["aldosterone"], color=COLORS["tertiary"], linewidth=1.5)
    ax11.set_ylabel("Normalized Level")
    ax11.set_title("Aldosterone")
    ax11.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    ax12 = fig.add_subplot(gs[3, 2])
    ax12.plot(t_days, ts["RSNA"], color=COLORS["quaternary"], linewidth=1.5)
    ax12.set_ylabel("Normalized Level")
    ax12.set_title("RSNA")
    ax12.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    # Row 5: ADH, ANP, Plasma Osm
    ax13 = fig.add_subplot(gs[4, 0])
    ax13.plot(t_days, ts["ADH"], color=COLORS["quinary"], linewidth=1.5)
    ax13.set_ylabel("Normalized Level")
    ax13.set_title("ADH (Vasopressin)")
    ax13.set_xlabel("Time (days)")
    ax13.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    ax14 = fig.add_subplot(gs[4, 1])
    ax14.plot(t_days, ts["ANP"], color=COLORS["senary"], linewidth=1.5)
    ax14.set_ylabel("Normalized Level")
    ax14.set_title("ANP")
    ax14.set_xlabel("Time (days)")
    ax14.axhline(1.0, color="gray", linestyle="--", alpha=0.5)

    ax15 = fig.add_subplot(gs[4, 2])
    ax15.plot(t_days, ts["plasma_osm"], color=COLORS["primary"], linewidth=1.5)
    ax15.set_ylabel("Osm (mOsm/kg)")
    ax15.set_title("Plasma Osmolality")
    ax15.set_xlabel("Time (days)")
    ax15.axhline(290, color="gray", linestyle="--", alpha=0.5)

    plt.savefig(os.path.join(output_dir, f"{safe_name}.png"))
    plt.close()
    print(f"  Plot saved: {safe_name}.png")


def run_and_plot_scenario(scenario_key, scenarios, output_dir):
    """Run a single scenario and generate plots."""
    sc = scenarios[scenario_key]
    print(f"\n▶ Running scenario: {sc['name']}...")

    params = sc["params"]
    t_span = sc["t_span"]

    # Time points for output (every 10 minutes)
    t_eval = np.arange(t_span[0], t_span[1], 10.0)

    try:
        sol = run_simulation(
            params=params,
            t_span=t_span,
            Na_intake_func=sc["Na_intake_func"],
            water_intake_func=sc["water_intake_func"],
            drug_effects_func=sc["drug_effects_func"],
            t_eval=t_eval,
            max_step=5.0,
        )

        if not sol.success:
            print(f"  ⚠ Solver warning: {sol.message}")

        # Extract time series
        ts = extract_timeseries(sol, params, sc["drug_effects_func"])

        # Print steady-state report
        print_steady_state_report(ts, sc["name"])

        # Generate plots
        plot_scenario(ts, sc, output_dir)

        return ts

    except Exception as e:
        print(f"  ✗ Error in scenario '{scenario_key}': {e}")
        import traceback
        traceback.print_exc()
        return None


def generate_comparison_plot(all_results, output_dir):
    """Generate a summary comparison across scenarios."""
    if not all_results:
        return

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("Cross-Scenario Comparison (Final Steady-State Values)",
                 fontsize=13, fontweight="bold")

    variables = [
        ("MAP", "MAP (mmHg)", 93),
        ("GFR", "GFR (mL/min)", 125),
        ("Na_excretion_daily", "Na Excretion (mEq/day)", 150),
        ("renin", "Renin (normalized)", 1.0),
        ("aldosterone", "Aldosterone (normalized)", 1.0),
        ("ADH", "ADH (normalized)", 1.0),
    ]

    names = list(all_results.keys())
    x = np.arange(len(names))

    for i, (var, label, normal) in enumerate(variables):
        ax = axes[i // 3, i % 3]
        values = [all_results[n][var][-1] if all_results[n] is not None else 0
                  for n in names]
        bars = ax.bar(x, values, color=list(COLORS.values())[i], alpha=0.8)
        ax.axhline(normal, color="red", linestyle="--", alpha=0.6, linewidth=1)
        ax.set_ylabel(label)
        ax.set_xticks(x)
        ax.set_xticklabels([n.replace("_", "\n") for n in names],
                           fontsize=7, rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "comparison_summary.png"))
    plt.close()
    print(f"\n  Comparison plot saved: comparison_summary.png")


def main():
    """Main entry point."""
    output_dir = "/home/claude/kidney_model/output"
    os.makedirs(output_dir, exist_ok=True)

    # Parse command line
    if "--list" in sys.argv:
        scenarios = get_all_scenarios()
        print("Available scenarios:")
        for key, sc in scenarios.items():
            print(f"  {key:20s} — {sc['name']}")
        return

    specific = None
    if "--scenario" in sys.argv:
        idx = sys.argv.index("--scenario")
        if idx + 1 < len(sys.argv):
            specific = sys.argv[idx + 1]

    # Initialize
    params = Parameters()
    scenarios = get_all_scenarios(params)

    if specific:
        if specific not in scenarios:
            print(f"Unknown scenario: {specific}")
            print(f"Available: {', '.join(scenarios.keys())}")
            return
        keys_to_run = [specific]
    else:
        keys_to_run = list(scenarios.keys())

    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   INTEGRATED KIDNEY MODEL — Simulation Engine              ║")
    print("║   6 Modules · ODE-based · Guyton/Karaaslan Framework       ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print(f"\nRunning {len(keys_to_run)} scenario(s)...")

    all_results = {}
    for key in keys_to_run:
        ts = run_and_plot_scenario(key, scenarios, output_dir)
        if ts is not None:
            all_results[key] = ts

    # Generate comparison
    if len(all_results) > 1:
        generate_comparison_plot(all_results, output_dir)

    print(f"\n✓ All outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
