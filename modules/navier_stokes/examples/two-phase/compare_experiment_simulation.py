# compare_experiment_simulation.py
"""
Compare experimental and simulation CSV data, generate plots, and (optionally) compute the Grid Convergence Index (GCI)
for mesh-refinement studies, with optional nondimensionalization of simulation data.

Quick start
-----------
Run the script with *no* flags to see the full help:

```bash
python compare_experiment_simulation.py --help
```

Typical usage
-------------
```bash
python compare_experiment_simulation.py \
    --exp       exp_front_height.csv \
    --sims      fine.csv medium.csv coarse.csv \
    --meshsizes 0.005 0.01 0.02 \
    --variable  compute_front_height \
    --output    comparison.png \
    [--tref 0.0005 --vref 0.06]
```

*   `--sims` should list files from **finest → coarsest** grid.
*   If ≥ 3 simulation grids are provided, the script will estimate the observed order of accuracy *p* and report GCI values per ASME V&V-20.
*   To nondimensionalize simulation time and the chosen variable, supply `--tref` (time scale) and `--vref` (variable scale). This does **not** alter experimental data.
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path
from typing import List, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ————————————————————————————————————————————————————————————————————————
# I/O & data preprocessing
# ————————————————————————————————————————————————————————————————————————

def read_csv(path: str | Path) -> pd.DataFrame:
    """Read a CSV file and trim whitespace/BOM from column names."""
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip().str.replace("\ufeff", "")
    return df


def align_on_time(
    dfs: List[pd.DataFrame], time_key: str = "time"
) -> List[pd.DataFrame]:
    """Interpolate DataFrames onto a common time vector."""
    time_union = np.unique(np.concatenate([df[time_key].to_numpy(dtype=float) for df in dfs]))
    aligned = []
    for df in dfs:
        interp = (
            df.set_index(time_key)
            .reindex(time_union)
            .interpolate("linear")
            .reset_index()
            .rename(columns={"index": time_key})
        )
        aligned.append(interp)
    return aligned


def nondimensionalize(
    sims: List[pd.DataFrame], time_key: str, var_key: str,
    tref: float, vref: float
) -> List[pd.DataFrame]:
    """Return new list of simulation DataFrames with nondimensional time & variable."""
    nd_sims: list[pd.DataFrame] = []
    for df in sims:
        df_nd = df.copy()
        df_nd[time_key] = df_nd[time_key] / tref
        df_nd[var_key]  = df_nd[var_key]  / vref
        nd_sims.append(df_nd)
    return nd_sims


# ————————————————————————————————————————————————————————————————————————
# GCI (Grid Convergence Index)
# ————————————————————————————————————————————————————————————————————————
_FS = 1.25  # ASME V&V-20 safety factor

def _richardson_p(
    phi1: np.ndarray, phi2: np.ndarray, phi3: np.ndarray,
    r21: float, r32: float
) -> float:
    eps = 1e-16
    num = np.log(np.linalg.norm(phi3 - phi2) / (np.linalg.norm(phi2 - phi1) + eps))
    den = np.log(r32 + eps)
    return num / den


def compute_gci(
    phi: List[np.ndarray], h: List[float], p: Optional[float] = None
) -> Dict[str, float]:
    if len(phi) < 3 or len(h) < 3:
        raise ValueError("Need at least three grids to compute GCI.")
    r21, r32 = h[1]/h[0], h[2]/h[1]
    if p is None:
        p = _richardson_p(phi[0], phi[1], phi[2], r21, r32)
    eps = 1e-16
    gci21 = _FS * np.linalg.norm(phi[1]-phi[0]) / (np.linalg.norm(phi[0])*(r21**p-1)+eps)
    gci32 = _FS * np.linalg.norm(phi[2]-phi[1]) / (np.linalg.norm(phi[1])*(r32**p-1)+eps)
    return {"p": p, "GCI21": gci21, "GCI32": gci32}

# ————————————————————————————————————————————————————————————————————————
# Plotting
# ————————————————————————————————————————————————————————————————————————

def plot_series_multi(
    exps: List[pd.DataFrame],
    exp_paths: List[str],
    sims: List[pd.DataFrame],
    sim_labels: List[str],
    variable: str,
    save_path: str | Path | None = None,
) -> None:
    """
+    Plot multiple experimental datasets (markers) and simulation curves (lines).
+    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # Plot each experiment as distinct markers
    for df, path in zip(exps, exp_paths):
        label = Path(path).stem
        ax.plot(df["time"], df[variable], marker='o', linestyle='',
                label=f"Exp: {label}")

    # Plot simulation curves
    for df, lbl in zip(sims, sim_labels):
        ax.plot(df["time"], df[variable], label=lbl)

    ax.set_xlabel("Time")
    ax.set_ylabel(variable)
    ax.grid(True, which="both", ls=":", lw=0.5)
    ax.legend()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    else:
        plt.show()

# ————————————————————————————————————————————————————————————————————————
# CLI
# ————————————————————————————————————————————————————————————————————————

def _build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="compare_experiment_simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """Compare experiment vs. simulation data,
            generate plots, compute GCI, and optionally nondimensionalize sims."""
        ),
        epilog="See --help for examples."
    )


def _parse_args(argv: List[str] | None=None) -> argparse.Namespace:
    p = _build_parser()
    p.add_argument(
        "--exps", nargs="+", required=True,
        help="Experiment CSV file paths (can specify multiple)"
    )
    p.add_argument(
        "--sims", nargs='+', required=True,
        help="Simulation CSV paths (finest → coarsest)"
    )
    p.add_argument(
        "--meshsizes", nargs='+', type=float, required=True,
        help="Mesh sizes for sims (fine → coarse)"
    )
    p.add_argument(
        "--variable", default="compute_front_height",
        help="Column name for quantity of interest"
    )
    p.add_argument(
        "--output", default=None,
        help="Optional output plot filename"
    )
    p.add_argument(
        "--tref", type=float, default=0.076326195,
        help="Reference time scale for nondimensionalization"
    )
    p.add_argument(
        "--vref", type=float, default=0.05715,
        help="Reference variable scale for nondimensionalization"
    )
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        p.print_help(sys.stderr)
        p.exit()
    return p.parse_args(argv)


def main(argv: List[str]|None=None) -> None:
    args = _parse_args(argv)
    if len(args.sims)!=len(args.meshsizes):
        raise SystemExit("--sims and --meshsizes length mismatch.")

    # Read experiment files
    exp_dfs = [read_csv(f) for f in args.exps]
    sim_dfs = [read_csv(f) for f in args.sims]

    # Validate variable column in each experiment
    var = args.variable
    for i, df in enumerate(exp_dfs):
        cols = df.columns.tolist()
        if var not in cols:
            alt = [c for c in cols if c.lower() == var.lower()]
            if alt:
                var = alt[0]
                print(f"Using '{var}' (case-insensitive) in experiment {args.exps[i]}")
            else:
                raise SystemExit(f"'{var}' not in experiment columns of {args.exps[i]}: {cols}")


    # Validate variable column in each simulation
    for i, df in enumerate(sim_dfs):
        if var not in df.columns.tolist():
            raise SystemExit(f"'{var}' not in simulation file {args.sims[i]}")

    # Optional nondimensionalization
    if args.tref and args.vref:
        sim_dfs = nondimensionalize(sim_dfs, "time", var, args.tref, args.vref)
        print(f"Nondimensionalized sims with tref={args.tref}, vref={args.vref}")

    # Align simulation data onto a common time vector
    aligned_sims = align_on_time(sim_dfs)

    # (Optional) align experiments as well if you want to compare on identical time points:
    # aligned_exps = align_on_time(exp_dfs)
    # But usually we just plot raw experimental markers.

    # Plot all experiments (markers) + simulations (lines)
    sim_labels = [f"Sim (h={h:g})" for h in args.meshsizes]
    plot_series_multi(exp_dfs, args.exps, aligned_sims, sim_labels, var, args.output)

    # Compute GCI if ≥3 meshes
    if len(aligned_sims) >= 3:
        idx = np.argsort(args.meshsizes)
        h_sorted = np.array(args.meshsizes)[idx]
        sims_sorted = [aligned_sims[i] for i in idx]
        phi = [df[var].to_numpy() for df in sims_sorted[:3]]
        gci = compute_gci(phi, h_sorted[:3])
        print("\n===== GCI =====")
        print(f"p = {gci['p']:.3f}, GCI21={gci['GCI21']:.3e}, GCI32={gci['GCI32']:.3e}")
    else:
        print("\n[GCI] Need 3+ sims; skipping.")


if __name__ == "__main__":
    main()
