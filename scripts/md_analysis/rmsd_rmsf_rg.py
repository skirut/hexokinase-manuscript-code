#!/usr/bin/env python3
"""RMSD, RMSF, and radius-of-gyration trajectory analysis.

Purpose
-------
Discovers Hexokinase simulation systems and calculates global RMSD,
per-residue RMSF, and radius of gyration from molecular-dynamics
trajectories. The script also performs the statistical processing and
publication-grid plotting implemented in the module.

Inputs
------
- Apo and/or bound systems under the configured data root.
- Compatible topology and trajectory files in each system directory.
- System names and pH values are parsed from the directory structure.

Outputs
-------
Under ``results/rmsd_rmsf_rg/``:
- analysis tables and intermediate data produced by the script
- publication-style RMSD/RMSF/Rg figures

Run from
--------
Run from the repository root with ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/md_analysis/rmsd_rmsf_rg.py``

Dependencies
------------
Python, NumPy, pandas, SciPy, Matplotlib, MDAnalysis, joblib, and
``scripts.utils.paths``.
"""

import os
import string
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from joblib import Parallel, delayed
from scipy.stats import ttest_rel, ttest_ind

plt.rcParams.update({
    "font.size": 9,
    "font.weight": "bold",          
    "axes.labelweight": "bold",     
    "axes.titleweight": "bold",     
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.titlesize": 12,
})

from scripts.utils.paths import APO_DIR, BOUND_DIR, RMSD_RMSF_RG_DIR

# External input/output paths
 
 
OUTDIR = RMSD_RMSF_RG_DIR
os.makedirs(OUTDIR, exist_ok=True)

ENZYME_NAMES = {
    "HK1": "Hexokinase 1",
    "HK2": "Hexokinase 2",
    "HK3": "Hexokinase 3",
    "HK4": "Hexokinase 4"
}

N_JOBS = 4
LIGAND_RESNAME = "G6P" 

CB_COLOR_APO = "#0072B2"  # Blue
CB_COLOR_BND = "#D55E00"  # Vermilion / Orange
CB_COLOR_LIG = "#009E73"  # Bluish Green

def ok_file(path):
    return path is not None and os.path.exists(path) and os.path.getsize(path) > 0

def load_universe(dirpath):
    if not dirpath:
        return None
    for subdir in ["prod_extended_500ns", "prod"]:
        top  = os.path.join(dirpath, subdir, "structure.prmtop")
        traj = os.path.join(dirpath, subdir, "output.xtc")
        if ok_file(top) and ok_file(traj):
            return mda.Universe(top, traj)
    return None

def parse_folder(name):
    parts = name.split('_')
    pdb = parts[0].lower()
    base = {"1v4s": "HK4", "2nzt": "HK2", "3hm8": "HK3", "1qha": "HK1"}.get(pdb, pdb.upper())
    
    ph = "7.0"
    for p in parts:
        if p.startswith("pH"):
            raw_ph = p.replace("pH", "")
            try:
                ph = f"{float(raw_ph):.1f}"
            except ValueError:
                ph = raw_ph
            
    mut = "WT"
    if len(parts) > 1:
        if parts[1] not in ["model", "G6P", f"pH{ph}"] and not parts[1].startswith("pH"):
            mut = parts[1]
            
    return base, mut, ph

def significance_stars(p):
    if p is None or np.isnan(p): return "ns"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"

def effect_size_cohens_d(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if x.size < 2 or y.size < 2: return np.nan
    sx2, sy2 = np.var(x, ddof=1), np.var(y, ddof=1)
    denom = np.sqrt((sx2 + sy2) / 2.0)
    if denom == 0: return np.nan
    return (np.mean(y) - np.mean(x)) / denom  

def analyze_one_system(apo_dir, bound_dir, base, mut, ph):
    warnings.filterwarnings("ignore")
    try:
        uA = load_universe(apo_dir)
        uB = load_universe(bound_dir)
        if uA is None and uB is None:
            return {"ok": False, "msg": f"missing top/traj in {base} {mut} pH{ph}"}

        sel_prot = "protein and name CA"
        
        # Apo calculation
        if uA is not None:
            protA = uA.select_atoms(sel_prot)
            RA = rms.RMSD(uA, select=sel_prot, ref_frame=0).run()
            tA_ns, pRMSD_A = RA.results.rmsd[:,1] / 1000.0, RA.results.rmsd[:,2]
            align.AlignTraj(uA, uA, select="protein and backbone", in_memory=True).run()
            RMSF_A = rms.RMSF(protA).run().results.rmsf
        else:
            tA_ns, pRMSD_A, RMSF_A = np.array([]), np.array([]), np.array([])

        # Bound calculation
        if uB is not None:
            protB = uB.select_atoms(sel_prot)
            RB = rms.RMSD(uB, select=sel_prot, ref_frame=0).run()
            tB_ns, pRMSD_B = RB.results.rmsd[:,1] / 1000.0, RB.results.rmsd[:,2]
            align.AlignTraj(uB, uB, select="protein and backbone", in_memory=True).run()
            RMSF_B = rms.RMSF(protB).run().results.rmsf

            ligB = uB.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
            if ligB.n_atoms > 0:
                RL = rms.RMSD(uB, select=f"resname {LIGAND_RESNAME} and not name H*", ref_frame=0).run()
                lRMSD, tl_ns = RL.results.rmsd[:,2], RL.results.rmsd[:,1] / 1000.0
            else:
                lRMSD, tl_ns = None, None
        else:
            tB_ns, pRMSD_B, RMSF_B = np.array([]), np.array([]), np.array([])
            lRMSD, tl_ns = None, None

        # Comparative statistics (only if both states exist)
        if RMSF_A.size > 0 and RMSF_B.size > 0 and RMSF_A.size == RMSF_B.size:
            dRMSF = RMSF_B - RMSF_A
        else:
            dRMSF = np.array([])

        n = min(len(pRMSD_A), len(pRMSD_B))
        if n >= 5:
            p_raw = ttest_rel(pRMSD_A[:n], pRMSD_B[:n]).pvalue if len(pRMSD_A) == len(pRMSD_B) else ttest_ind(pRMSD_A[:n], pRMSD_B[:n], equal_var=False).pvalue
            d_val = effect_size_cohens_d(pRMSD_A[:n], pRMSD_B[:n])
            stars_val = significance_stars(p_raw)
        else:
            p_raw, d_val, stars_val = np.nan, np.nan, "N/A"

        return {
            "ok": True, "base": base, "mut": mut, "pH": ph,
            "time_ns_apo": tA_ns, "prot_rmsd_apo": pRMSD_A,
            "time_ns_bnd": tB_ns, "prot_rmsd_bnd": pRMSD_B,
            "lig_time_ns": tl_ns, "lig_rmsd": lRMSD,
            "rmsf_apo": RMSF_A, "rmsf_bnd": RMSF_B, "delta_rmsf": dRMSF,
            "p_raw": p_raw, "cohens_d": d_val, "stars": stars_val
        }
    except Exception as e:
        return {"ok": False, "base": base, "mut": mut, "pH": ph, "msg": str(e)}

def discover_systems():
    apo_dict = {}
    if os.path.exists(APO_DIR):
        for name in sorted(os.listdir(APO_DIR)):
            p = os.path.join(APO_DIR, name)
            if os.path.isdir(p):
                parsed = parse_folder(name)
                apo_dict[parsed] = p

    bound_dict = {}
    if os.path.exists(BOUND_DIR):
        for name in sorted(os.listdir(BOUND_DIR)):
            p = os.path.join(BOUND_DIR, name)
            if os.path.isdir(p):
                parsed = parse_folder(name)
                bound_dict[parsed] = p

    all_keys = sorted(list(set(apo_dict.keys()) | set(bound_dict.keys())))
    systems = []
    for key in all_keys:
        apo_path = apo_dict.get(key, None)
        bnd_path = bound_dict.get(key, None)
        systems.append((apo_path, bnd_path, key[0], key[1], key[2]))
    return systems

def plot_protein_grid(df_prot, base, outprefix):
    df = df_prot.copy()
    df["pH_float"] = df["pH"].astype(float)
    
    unique_phs = sorted(df["pH_float"].unique())
    muts = df["mut"].unique().tolist()
    if "WT" in muts:
        muts.remove("WT")
        mut_order = ["WT"] + sorted(muts)
    else:
        mut_order = sorted(muts)

    n_muts = len(mut_order)
    n_phs = len(unique_phs)
    is_2x2_mode = (n_muts == 1 and n_phs <= 4)
    
    if is_2x2_mode:
        nrows, ncols_pairs = 2, 2
        fig, axes = plt.subplots(nrows, ncols_pairs * 2, figsize=(10, 5.5))
        fig.subplots_adjust(left=0.10, right=0.98, top=0.78, bottom=0.15, wspace=0.45, hspace=0.85)
    else:
        nrows, ncols_pairs = n_muts, n_phs
        fig, axes = plt.subplots(nrows, ncols_pairs * 2, figsize=(5.0 * n_phs, 2.4 * nrows))
        axes = np.atleast_2d(axes)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.12, wspace=0.45, hspace=0.65)

    letters = string.ascii_uppercase
    enz_name = ENZYME_NAMES.get(base, base)
    fig.suptitle(f"RMSD & RMSF Analysis for {base.upper()} ({enz_name})", fontsize=18, fontweight="bold", y=0.96)
    
    for r in range(nrows):
        for c in range(ncols_pairs):
            ax_rmsd = axes[r, c*2]
            ax_rmsf = axes[r, c*2 + 1]
            
            if is_2x2_mode:
                mut = mut_order[0]
                ph_idx = r * 2 + c
                if ph_idx >= n_phs:
                    ax_rmsd.axis('off'); ax_rmsf.axis('off')
                    continue
                ph = unique_phs[ph_idx]
                label_letter = letters[ph_idx]
            else:
                mut = mut_order[r]
                ph = unique_phs[c]
                label_letter = letters[r]
                
            subset = df[(df["mut"] == mut) & (df["pH_float"] == ph)]
            
            if is_2x2_mode:
                ax_rmsd.text(-0.30, 1.15, f"{label_letter})", transform=ax_rmsd.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')
                ax_rmsd.text(1.0, 1.15, f"pH {ph:.1f}", transform=ax_rmsd.transAxes, fontsize=13, fontweight="bold", ha="center", va="bottom")
            else:
                if r == 0:
                    ax_rmsd.text(1.0, 1.15, f"pH {ph:.1f}", transform=ax_rmsd.transAxes, fontsize=13, fontweight="bold", ha="center", va="bottom")
                if c == 0:
                    ax_rmsd.text(-0.30, 0.5, label_letter if r < len(letters) else "", transform=ax_rmsd.transAxes, fontsize=16, fontweight='bold', va='center', ha='right')

            if subset.empty:
                ax_rmsd.axis('off'); ax_rmsf.axis('off')
                continue
                
            row = subset.iloc[0]
            title_str = f"RMSD of {base} ({enz_name}) [{row['mut']}]"

            ta, ya = np.asarray(row["time_ns_apo"]), np.asarray(row["prot_rmsd_apo"])
            tb, yb = np.asarray(row["time_ns_bnd"]), np.asarray(row["prot_rmsd_bnd"])
            
            if ta.size > 0: ax_rmsd.plot(ta, ya, color=CB_COLOR_APO, linewidth=1.2, label="Protein (apo)")
            if tb.size > 0: ax_rmsd.plot(tb, yb, color=CB_COLOR_BND, linewidth=1.2, label="Protein (bound)")
            
            if row["lig_rmsd"] is not None:
                tl, yl = np.asarray(row["lig_time_ns"]), np.asarray(row["lig_rmsd"])
                if tl.size > 0: ax_rmsd.plot(tl, yl, linestyle="--", linewidth=1.2, color=CB_COLOR_LIG, label=f"Ligand ({LIGAND_RESNAME})")
            
            ax_rmsd.set_title(title_str, fontsize=9, fontweight="bold")
            ax_rmsd.set_ylabel("RMSD (nm)", fontweight="bold")
            ax_rmsd.set_xlabel("Time (ns)", fontweight="bold")
            ax_rmsd.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
            ax_rmsd.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
            ax_rmsd.grid(True, alpha=0.3)
            
            if row["stars"] != "N/A" and not np.isnan(row["p_raw"]):
                stat_txt = f"{row['stars']}  p={row['p_raw']:.2e}"
                if not np.isnan(row["cohens_d"]): stat_txt += f", d={row['cohens_d']:+.2f}"
            else:
                stat_txt = "Apo State Only"
            
            ax_rmsd.text(0.97, 0.05, stat_txt, transform=ax_rmsd.transAxes,
                    ha="right", va="bottom", fontsize=7, fontweight="bold", zorder=10,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=1.0))

            ra, rb, dr = np.asarray(row["rmsf_apo"]), np.asarray(row["rmsf_bnd"]), np.asarray(row["delta_rmsf"])
            if ra.size > 0: ax_rmsf.plot(ra, color=CB_COLOR_APO, linewidth=1.2, label="Protein RMSF (apo)")
            if rb.size > 0: ax_rmsf.plot(rb, color=CB_COLOR_BND, linewidth=1.2, label="Protein RMSF (bound)")
            if dr.size > 0: ax_rmsf.plot(dr, color=CB_COLOR_LIG, linewidth=1.2, linestyle="--", label="ΔRMSF (bound−apo)")
            
            ax_rmsf.set_title("RMSF Comparison", fontsize=9, fontweight="bold")
            ax_rmsf.set_ylabel("RMSF (nm)", fontweight="bold")
            ax_rmsf.set_xlabel("Residue index", fontweight="bold")
            ax_rmsf.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
            ax_rmsf.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
            ax_rmsf.grid(True, alpha=0.3)

            if (is_2x2_mode and ph_idx == 0) or (not is_2x2_mode and r == 0 and c == 0):
                current_rmsd_top = ax_rmsd.get_ylim()[1]
                current_rmsf_top = ax_rmsf.get_ylim()[1]
                ax_rmsd.set_ylim(top=current_rmsd_top * 1.4)
                ax_rmsf.set_ylim(top=current_rmsf_top * 1.4)
                
                leg_rmsd = ax_rmsd.legend(fontsize=7, framealpha=1.0, loc="upper left", edgecolor="none")
                if leg_rmsd: leg_rmsd.set_zorder(10)
                leg_rmsf = ax_rmsf.legend(fontsize=7, framealpha=1.0, loc="upper left", edgecolor="none")
                if leg_rmsf: leg_rmsf.set_zorder(10)
    
    legend_text = (
        "P-value stars: ns (p > 0.05), * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001) | "
        "Cohen's d: ~0.2 small, ~0.5 medium, ~0.8 large | sign: + (bound > apo), − (apo > bound)"
    )
    y_pos = 0.01 if is_2x2_mode else 0.02
    fig.text(0.5, y_pos, legend_text, ha="center", va="center", fontsize=9, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none"))

    out_tiff, out_png = f"{outprefix}.tiff", f"{outprefix}.png"
    fig.savefig(out_tiff, dpi=300, facecolor="white", format="tiff", bbox_inches="tight")
    fig.savefig(out_png, dpi=300, facecolor="white", format="png", bbox_inches="tight")
    plt.close(fig)
    print(f"[saved publication grid for {base}] -> {out_tiff}, {out_png}")

def main():
    systems = discover_systems()
    print(f"[INFO] Found {len(systems)} total hexokinase systems")
    
    if not systems:
        print("[INFO] nothing to do"); return

    results = Parallel(n_jobs=N_JOBS)(
        delayed(analyze_one_system)(apo, bnd, base, mut, ph) for (apo, bnd, base, mut, ph) in systems
    )

    rows = [r for r in results if r.get("ok", False)]
    if not rows:
        print("[warn] No successful analyses; exiting.")
        return

    df = pd.DataFrame(rows)
    df.to_pickle(os.path.join(OUTDIR, "rmsd_rmsf_pairs.pkl"))
    
    df_stats = df[["base","mut","pH","p_raw","cohens_d","stars"]].copy()
    df_stats.to_csv(os.path.join(OUTDIR, "rmsd_rmsf_stats.csv"), index=False)

    for base, sub in df.groupby("base"):
        outprefix = os.path.join(OUTDIR, f"MASTER_RMSD_RMSF_GRID_PUBLICATION_{base.upper()}")
        plot_protein_grid(sub, base, outprefix)

if __name__ == "__main__":
    main()

