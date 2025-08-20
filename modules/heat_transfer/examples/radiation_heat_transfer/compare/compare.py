import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple

def infer_cols(df: pd.DataFrame) -> Tuple[str, str]:
    lower = {c.lower(): c for c in df.columns}
    x_candidates = ["x", "r", "position", "pos", "radius"]
    t_candidates = ["t", "temperature", "temp"]
    x_col = next((lower[c] for c in x_candidates if c in lower), None)
    t_col = next((lower[c] for c in t_candidates if c in lower), None)
    if x_col is None or t_col is None:
        num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(num) < 2:
            raise ValueError("Cannot infer x/T columns and <2 numeric columns present.")
        x_col, t_col = num[0], num[1]
    return x_col, t_col

def prep_ref(df_ref: pd.DataFrame, x_col: str, t_col: str) -> pd.DataFrame:
    ref = df_ref[[x_col, t_col]].copy()
    ref[x_col] = pd.to_numeric(ref[x_col], errors="coerce")
    ref[t_col] = pd.to_numeric(ref[t_col], errors="coerce")
    ref = ref.dropna().sort_values(x_col)
    ref = ref.groupby(x_col, as_index=False)[t_col].mean()
    return ref.reset_index(drop=True)

def interp_bracketed(x_ref: np.ndarray, t_ref: np.ndarray, xq: float) -> float:
    n = x_ref.size
    if n < 2:
        raise ValueError("Reference needs at least two distinct x values.")
    i = np.searchsorted(x_ref, xq, side="left")
    if i <= 0:
        i0, i1 = 0, 1
    elif i >= n:
        i0, i1 = n - 2, n - 1
    else:
        i0, i1 = i - 1, i
    x0, x1 = x_ref[i0], x_ref[i1]
    t0, t1 = t_ref[i0], t_ref[i1]
    if np.isclose(x0, x1):
        return float(t0)
    w = (xq - x0) / (x1 - x0)
    return float(t0 + w * (t1 - t0))

def main(ref_path: str, moose_path: str, out_dir: str = ".", plot_path: str = None):
    df_ref_raw = pd.read_csv(ref_path)
    df_moose_raw = pd.read_csv(moose_path)
    ref_x, ref_t = infer_cols(df_ref_raw)
    moose_x, moose_t = infer_cols(df_moose_raw)
    df_ref = prep_ref(df_ref_raw, ref_x, ref_t)
    x_ref = df_ref[ref_x].to_numpy(float)
    t_ref = df_ref[ref_t].to_numpy(float)
    x_moose = pd.to_numeric(df_moose_raw[moose_x], errors="coerce").to_numpy(float)
    t_moose = pd.to_numeric(df_moose_raw[moose_t], errors="coerce").to_numpy(float)
    mask_valid = np.isfinite(x_moose) & np.isfinite(t_moose)
    x_moose = x_moose[mask_valid]
    t_moose = t_moose[mask_valid]
    t_ref_est = np.array([interp_bracketed(x_ref, t_ref, xv) for xv in x_moose], dtype=float)
    err_abs = t_moose - t_ref_est
    eps = 1e-12
    err_rel_pct = np.where(np.abs(t_ref_est) > eps, err_abs / t_ref_est * 100.0, np.nan)
    ref_base = os.path.splitext(os.path.basename(ref_path))[0]
    moose_base = os.path.splitext(os.path.basename(moose_path))[0]
    out_path = os.path.join(out_dir, f"compare_results/compare_{ref_base}_{moose_base}.csv")
    out = pd.DataFrame({
        "x": x_moose,
        "T_moose": t_moose,
        "T_ref_interp": t_ref_est,
        "err_abs": err_abs,
        "err_rel_pct": err_rel_pct
    })
    out.to_csv(out_path, index=False)
    if plot_path is None:
        plot_path = os.path.join(out_dir, f"compare_results/error_vs_x_{moose_base}.png")
    plt.figure()
    plt.plot(x_moose, err_abs, linestyle="none", marker="o", markersize=6)
    plt.xlabel("x")
    plt.ylabel("Absolute error (T_moose - T_ref_interp)")
    plt.title("Error vs. position")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    with np.errstate(invalid="ignore"):
        idx_rel = np.nanargmax(np.abs(err_rel_pct)) if np.any(~np.isnan(err_rel_pct)) else None
    idx_abs = int(np.argmax(np.abs(err_abs))) if err_abs.size > 0 else None
    if idx_rel is not None:
        print(f"\nMax relative error: {err_rel_pct[idx_rel]:.6g} %")
        print(f"At x = {x_moose[idx_rel]:.6g}")
        print(f"T_moose={t_moose[idx_rel]:.6g}, T_ref={t_ref_est[idx_rel]:.6g}, abs_err={err_abs[idx_rel]:.6g}")
    else:
        print("\nMax relative error: N/A (reference near zero)")
    if idx_abs is not None:
        print(f"\nMax absolute error: {np.abs(err_abs[idx_abs]):.6g}")
        print(f"At x = {x_moose[idx_abs]:.6g}")
        print(f"T_moose={t_moose[idx_abs]:.6g}, T_ref={t_ref_est[idx_abs]:.6g}, rel_err_pct={err_rel_pct[idx_abs]:.6g}")
    print(f"\nCSV : {out_path}")
    print(f"Plot : {plot_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
            print("Usage: python your_script_name.py <file1> <file2>")
    else:
        main(sys.argv[1], sys.argv[2])

