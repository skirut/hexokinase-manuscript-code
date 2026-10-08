#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Integrated G6P hydrogen-bond, center-of-mass, and contact analysis.

Purpose
-------
Analyzes protein-G6P hydrogen bonds, ligand-protein center-of-mass
distance, and binding-site contacts for the bound Hexokinase systems.
The script also builds summary tables and publication-style time-series
and density plots across systems and pH conditions.

Inputs
------
- Bound systems under ``$HEXOKINASE_DATA_ROOT/G6P_structure``
- ``structure.prmtop`` and ``output.xtc`` (or supported trajectory
  alternatives) within the system production directories.
- G6P ligand residue identification is handled by the script.

Outputs
-------
Under ``results/ligand_interactions/``:
- ``hbond/`` hydrogen-bond outputs
- ``com/`` center-of-mass outputs
- ``contacts/`` binding-site contact outputs
- summary tables and publication figures

Run from
--------
Run from the repository root with ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/ligand_interactions/hbond_com_contacts.py``

Dependencies
------------
Python, NumPy, pandas, SciPy, Matplotlib, and ``scripts.utils.paths``.
"""


import os
import re
import string
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scripts.utils.paths import BOUND_DIR as BOUND_DATA_DIR, LIGAND_INTERACTIONS_DIR
import matplotlib.ticker as ticker
from scipy.stats import gaussian_kde

# ==========================================
# 1. BOLD MATPLOTLIB GLOBAL CONFIGURATION
# ==========================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "font.weight": "bold",          
    "axes.labelweight": "bold",     
    "axes.titleweight": "bold",     
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
})

# ===================== USER SETTINGS =====================
BOUND_DIR = BOUND_DATA_DIR
OUTROOT = LIGAND_INTERACTIONS_DIR
OUTROOT_HB = os.path.join(OUTROOT, "hbond")
OUTROOT_CO = os.path.join(OUTROOT, "com")
OUTROOT_CT = os.path.join(OUTROOT, "contacts")

TOPO_NAME  = "structure.prmtop"
TRAJ_NAME  = "output.xtc"

STRIDE        = 1

HB_DIST_MAIN  = 3.5
HB_ANGLE_MAIN = 150.0
USE_PERIODIC  = True
LIG_PREFERRED = "G6P"

CONSERVATIVE          = False      
HB_DIST_CONSERVATIVE  = 3.2
DONORS_N_ONLY         = True
COUNT_MODE = "atom"   

BS_DIST_A        = 5.0   
CONTACT_CUTOFF_A = 4.5   

MAX_WORKERS = min(20, os.cpu_count() or 4)

CBLUE   = "#0072B2"
CORANGE = "#E69F00"
CGREEN  = "#009E73"
CRED    = "#D55E00"
CPURPLE = "#CC79A7"
CGRAY   = "#999999"

ENZYME_NAME = {"1qha": "HK1", "2nzt": "HK2", "3hm8": "HK3"}

ROLL_WIN = 50

# System mappings based on folder names
SYSTEM_MAPPING = {
    '1qha': '1QHA_model_with_G6P',
    '1qha_G87E': '1QHA_G87E_G6P',
    '2nzt': '2nzt_model_with_G6P',
    '2nzt_D84Y': '2nzt_D84Y_G6P',
    '2nzt_D861Y': '2nzt_D861Y_G6P',
    '2nzt_T680M': '2nzt_T680M_G6P',
    '3hm8': '3hm8_model_with_G6P',
    '3hm8_D867N': '3hm8_D867N_G6P',
    '3hm8_D867Y': '3hm8_D867Y_G6P',
    '3hm8_T546M': '3hm8_T546M_G6P'
}
PH_LEVELS = ['3.0', '5.0', '7.0']

# ========================================================

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

ION_NAMES   = {"NA","K","CL","MG","CA","ZN","MN","CU","FE","CO","F","BR","IOD","CS","SR","LI","RB","SOD","POT"}
WATER_NAMES = {"HOH","WAT","TIP3","SOL","SPC","T3P","T4P"}

def ensure_outdirs():
    for root in (OUTROOT_HB, OUTROOT_CO, OUTROOT_CT):
        os.makedirs(root, exist_ok=True)
        os.makedirs(os.path.join(root, "per_system"), exist_ok=True)
        os.makedirs(os.path.join(root, "master_grids"), exist_ok=True)

def list_bound_systems():
    systems = []
    for internal_name, bound_prefix in SYSTEM_MAPPING.items():
        parts = internal_name.split('_')
        base = parts[0].lower()
        mut = parts[1] if len(parts) > 1 else 'WT'
        rowkey = f"{base}_{mut}"

        for ph in PH_LEVELS:
            ph_float = float(ph)
            sysdir = os.path.join(BOUND_DIR, f"{bound_prefix}_pH{ph}", "prod")
            topo = os.path.join(sysdir, TOPO_NAME)
            traj = os.path.join(sysdir, TRAJ_NAME)

            if os.path.isfile(topo) and os.path.isfile(traj):
                systems.append((rowkey, ph_float, sysdir, topo, traj))
    return systems

def detect_ligand_resname(u, preferred=LIG_PREFERRED):
    try:
        if u.select_atoms(f"resname {preferred}").n_atoms > 0: return preferred
    except Exception: pass
    for res in u.residues:
        rn = (res.resname or "").upper()
        if rn in WATER_NAMES or rn in ION_NAMES: continue
        if res.atoms.select_atoms("protein").n_atoms > 0: continue
        if 2 <= res.atoms.n_atoms <= 400: return rn
    return None

def donors_acceptors(prot_sel, lig_sel, conservative=False):
    if conservative:
        return (prot_sel.select_atoms("name N* or name O*"), prot_sel.select_atoms("name N* or name O*"),
                lig_sel.select_atoms("name N* or name O*"), lig_sel.select_atoms("name N* or name O*"), HB_DIST_CONSERVATIVE)
    return (prot_sel.select_atoms("name N* or name O* or name S*"), prot_sel.select_atoms("name N* or name O*"),
            lig_sel.select_atoms("name N* or name O* or name S*"), lig_sel.select_atoms("name N* or name O*"), HB_DIST_MAIN)

def run_hbond_analysis(u, lig_resname, stride, d_cut, angle_cut):
    cols = ["frame","donor_idx","hydrogen_idx","acceptor_idx","distance","angle"]
    try:
        from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
        h = HBA(universe=u, between=["protein", f"resname {lig_resname}"], 
                d_a_cutoff=d_cut, d_h_a_angle_cutoff=angle_cut, update_selections=False)
        h.run(step=stride)
        arr = getattr(h.results, "hbonds", None) if hasattr(h, "results") else None
        return pd.DataFrame(arr, columns=cols) if arr is not None else pd.DataFrame(columns=cols)
    except Exception:
        pass
    try:
        from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis as HBA_legacy
        h = HBA_legacy(u, "protein", f"resname {lig_resname}", d_cut, angle_cut)
        h.run(step=stride)
        arr = getattr(getattr(h, "results", h), "hbonds", None)
        if arr is None: return pd.DataFrame(columns=cols)
        df = pd.DataFrame(arr, columns=cols[:arr.shape[1]])
        for c in cols[df.shape[1]:]: df[c] = np.nan
        return df[cols]
    except Exception:
        return pd.DataFrame(columns=cols)

def hbond_proxy_counts(u, D_prot, A_prot, D_lig, A_lig, d_cut, stride=1, mode="atom"):
    from MDAnalysis.lib.distances import distance_array
    analyzed = list(range(0, len(u.trajectory), stride))
    counts, hits = np.zeros(len(analyzed), dtype=int), []
    for k, fr in enumerate(analyzed):
        u.trajectory[fr]
        c, frame_pairs = 0, set()
        for D, A in [(D_prot, A_lig), (D_lig, A_prot)]:
            if D.n_atoms and A.n_atoms:
                di, ai = (distance_array(D.positions, A.positions) <= d_cut).nonzero()
                if mode == "atom":
                    c += len(di)
                    for d, a in zip(di, ai): hits.append((f"{D.atoms[d].resname}:{D.atoms[d].name}→{A.atoms[a].name}", k))
                else:
                    for d, a in zip(di, ai): frame_pairs.add((D.atoms[d].resid, A.atoms[a].resid))
        if mode == "residue":
            c = len(frame_pairs)
            for dp, ap in frame_pairs: hits.append((f"Res{dp}→Res{ap}", k))
        counts[k] = c
    start = int(0.8 * len(analyzed)) if len(analyzed) >= 10 else 0
    top_pairs = [(p, cnt / max(1, len(analyzed) - start) * 100.0) for p, cnt in Counter([p for p, k in hits if k >= start]).most_common(10)]
    return counts, top_pairs

def compute_binding_site(u, lig_ag, bs_dist_A=BS_DIST_A):
    u.trajectory[0]
    try: return u.select_atoms(f"protein and around {bs_dist_A} resname {lig_ag.residues.resnames[0]}")
    except Exception: return u.select_atoms("protein")

def compute_com_and_contacts(u, lig_ag, site_ag, contact_cutoff_A=CONTACT_CUTOFF_A, stride=1):
    from MDAnalysis.lib.distances import distance_array
    times_ps, com_nm, contacts = [], [], []
    lig_heavy = lig_ag.select_atoms("not name H*") if lig_ag.n_atoms else lig_ag
    site_heavy = site_ag.select_atoms("not name H*") if site_ag.n_atoms else site_ag
    for ts in u.trajectory[::stride]:
        times_ps.append(float(getattr(ts, "time", ts.frame * (getattr(ts, "dt", 1.0) or 1.0))))
        try: com_nm.append(np.linalg.norm(site_heavy.center_of_mass() - lig_heavy.center_of_mass()) / 10.0)
        except Exception: com_nm.append(np.nan)
        contacts.append(int((distance_array(lig_heavy.positions, site_heavy.positions) <= contact_cutoff_A).sum()) if site_heavy.n_atoms else 0)
    return np.asarray(times_ps) / 1000.0, np.asarray(com_nm), np.asarray(contacts, dtype=int)

def smooth_series(vals, window=ROLL_WIN):
    return pd.Series(vals).rolling(window=max(3, min(window, len(vals))), min_periods=1).mean().to_numpy() if vals.size else vals

def analyze_one(args):
    rowkey, pH, sysdir, topo, traj = args
    import MDAnalysis as mda
    out_hb = os.path.join(OUTROOT_HB, "per_system", f"{rowkey}_pH{pH:g}")
    out_co = os.path.join(OUTROOT_CO, "per_system", f"{rowkey}_pH{pH:g}")
    out_ct = os.path.join(OUTROOT_CT, "per_system", f"{rowkey}_pH{pH:g}")
    for d in (out_hb, out_co, out_ct): os.makedirs(d, exist_ok=True)
    try: u = mda.Universe(topo, traj)
    except Exception as e: return {"ok": False, "rowkey": rowkey, "pH": pH, "msg": str(e)}

    lig_resname = detect_ligand_resname(u, LIG_PREFERRED)
    if not lig_resname: return {"ok": False, "rowkey": rowkey, "pH": pH, "msg": "no ligand"}

    prot = u.select_atoms("protein")
    lig  = u.select_atoms(f"resname {lig_resname}")

    # H-Bonds
    conservative = bool(CONSERVATIVE)
    D_prot, A_prot, D_lig, A_lig, d_cut = donors_acceptors(prot, lig, conservative=conservative)
    hb_df = run_hbond_analysis(u, lig_resname, STRIDE, d_cut, HB_ANGLE_MAIN)
    method = "HBA" if hb_df is not None and not hb_df.empty else "heavy-only"
    
    t_ns = np.asarray([float(getattr(ts, "time", ts.frame * (getattr(ts, "dt", 1.0) or 1.0))) for ts in u.trajectory[::STRIDE]]) / 1000.0
    if method == "HBA":
        counts_full = np.zeros(int(hb_df["frame"].max()) + 1, dtype=int) if not hb_df.empty else []
        for fr in hb_df["frame"].astype(int):
            if 0 <= fr < len(counts_full): counts_full[fr] += 1
        counts = np.array([counts_full[i * STRIDE] if (i * STRIDE) < len(counts_full) else 0 for i in range(len(t_ns))])
    else:
        counts, _ = hbond_proxy_counts(u, D_prot, A_prot, D_lig, A_lig, d_cut, stride=STRIDE, mode=COUNT_MODE)
    
    hb_ma = smooth_series(counts)
    hb_csv = os.path.join(out_hb, "hbonds_counts.csv")
    pd.DataFrame({"time_ns": t_ns, "count": counts, "count_ma": hb_ma}).to_csv(hb_csv, index=False)

    # COM / Contacts
    site = compute_binding_site(u, lig, bs_dist_A=BS_DIST_A)
    _, com_nm, contacts = compute_com_and_contacts(u, lig, site, contact_cutoff_A=CONTACT_CUTOFF_A, stride=STRIDE)
    
    com_csv = os.path.join(out_co, "com_counts.csv")
    pd.DataFrame({"time_ns": t_ns, "com_nm": com_nm, "com_nm_ma": smooth_series(com_nm)}).to_csv(com_csv, index=False)
    
    ct_csv = os.path.join(out_ct, "contacts_counts.csv")
    pd.DataFrame({"time_ns": t_ns, "contacts": contacts, "contacts_ma": smooth_series(contacts)}).to_csv(ct_csv, index=False)

    return {
        "ok": True, "rowkey": rowkey, "pH": pH, "ligand": lig_resname,
        "hb_csv": hb_csv, "com_csv": com_csv, "ct_csv": ct_csv
    }

def build_summary(results):
    out = []
    for r in results:
        if not r.get("ok"): continue
        parts = r["rowkey"].split("_")
        base_root = parts[0].lower()
        mut = parts[1] if len(parts) > 1 else "WT"
        out.append({
            "rowkey": r["rowkey"], "base": base_root, "mut": mut, 
            "pH": r["pH"], "hb_csv": r["hb_csv"], "com_csv": r["com_csv"], "ct_csv": r["ct_csv"]
        })
    return pd.DataFrame(out).sort_values(["base", "mut", "pH"])

# =========================
# PUBLICATION GRIDS 
# =========================
def plot_publication_grid(df_prot, base, metric):
    unique_phs = sorted(df_prot["pH"].unique())
    muts = df_prot["mut"].unique().tolist()
    if "WT" in muts:
        muts.remove("WT")
        mut_order = ["WT"] + sorted(muts)
    else: mut_order = sorted(muts)

    n_muts, n_phs = len(mut_order), len(unique_phs)
    
    nrows, ncols_pairs = n_muts, n_phs
    fig, axes = plt.subplots(nrows, ncols_pairs * 2, figsize=(3.5 * ncols_pairs * 2, 3.2 * nrows))
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.12, wspace=0.35, hspace=0.65)

    letters = string.ascii_uppercase
    enz_name = ENZYME_NAME.get(base, base.upper())
    
    # Configuration based on metric
    if metric == "hbond":
        y_col, ma_col, csv_col = "count", "count_ma", "hb_csv"
        c_raw, c_ma = CBLUE, CORANGE
        ylabel, title_prefix = "H-Bonds", "HB"
        metric_title = "Hydrogen Bonds"
        banner_text = "Hydrogen Bonds (Count) | Values inside black boxes: Mean ± SD"
    elif metric == "com":
        y_col, ma_col, csv_col = "com_nm", "com_nm_ma", "com_csv"
        c_raw, c_ma = CGREEN, CPURPLE
        ylabel, title_prefix = "COM Distance (nm)", "COM"
        metric_title = "Center of Mass (COM)"
        banner_text = "Ligand-Site COM Distance (nm) | Lower values → Closer | Values inside black boxes: Mean ± SD"
    else: 
        y_col, ma_col, csv_col = "contacts", "contacts_ma", "ct_csv"
        c_raw, c_ma = CRED, CBLUE
        ylabel, title_prefix = "Contacts", "Contacts"
        metric_title = "Heavy-Atom Contacts"
        banner_text = "Ligand Heavy-Atom Contacts (Count) | Higher values → More interactions | Values inside black boxes: Mean ± SD"

    # MAIN FIGURE TITLE INJECTION
    fig.suptitle(f"{metric_title} Analysis for {base.upper()} ({enz_name})", fontsize=18, fontweight="bold", y=0.96)

    for r in range(nrows):
        for c in range(ncols_pairs):
            ax_ts, ax_k = axes[r, c*2], axes[r, c*2 + 1]
            mut, ph, label_letter = mut_order[r], unique_phs[c], letters[r]
            
            if r == 0: 
                ax_ts.text(1.05, 1.25, f"pH {ph:.1f}", transform=ax_ts.transAxes, fontsize=16, fontweight="bold", ha="center", va="bottom")
            if c == 0: 
                ax_ts.text(-0.35, 0.5, label_letter if r < len(letters) else "", transform=ax_ts.transAxes, fontsize=20, fontweight='bold', va='center', ha='right')

            subset = df_prot[(df_prot["mut"] == mut) & (df_prot["pH"] == ph)]
            if subset.empty:
                ax_ts.axis('off'); ax_k.axis('off'); continue
                
            df_csv = pd.read_csv(subset.iloc[0][csv_col])
            t_ns, raw_val, ma_val = df_csv["time_ns"].values, df_csv[y_col].values, df_csv[ma_col].values
            
            # Filter valid raw values, dropping NaNs to prevent statistical crashes
            valid_raw = raw_val[~np.isnan(raw_val)]

            # Time Series (using raw_val so NaNs correctly render as gaps on the plot)
            ax_ts.plot(t_ns, raw_val, color=c_raw, alpha=0.4, linewidth=1.0, label="Raw")
            ax_ts.plot(t_ns, ma_val, color=c_ma, linewidth=1.5, label="Running Avg")
            
            ax_ts.set_title(f"{title_prefix} of {enz_name} [{mut}]", fontsize=10, fontweight="bold")
            ax_ts.set_ylabel(ylabel, fontweight="bold")
            ax_ts.set_xlabel("Time (ns)", fontweight="bold")
            ax_ts.set_xlim(0, 300)
            ax_ts.xaxis.set_major_locator(ticker.MultipleLocator(50))
            ax_ts.grid(True, alpha=0.3)
            
            # Escape route for entirely NaN data
            if len(valid_raw) == 0:
                ax_ts.text(0.5, 0.5, "No Valid Data\n(Selection Empty)", transform=ax_ts.transAxes, ha="center", va="center", fontweight="bold", color="red")
                ax_k.text(0.5, 0.5, "No Valid Data", transform=ax_k.transAxes, ha="center", va="center", fontweight="bold", color="red")
                continue

            avg = float(np.mean(valid_raw))
            std = float(np.std(valid_raw, ddof=1) if len(valid_raw) > 1 else 0.0)
            ax_ts.text(0.97, 0.05, f"{avg:.1f} ± {std:.1f}", transform=ax_ts.transAxes, ha="right", va="bottom", fontsize=8, fontweight="bold", zorder=10, bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1.5, alpha=1.0))

            # Density Plot (Using only valid_raw)
            try:
                if len(np.unique(valid_raw)) > 1:
                    density = gaussian_kde(valid_raw)
                    xs = np.linspace(np.min(valid_raw), np.max(valid_raw), 200)
                    ax_k.plot(xs, density(xs), label="Density (KDE)", color=c_raw, linewidth=1.5)
                else:
                    raise ValueError("Not enough variance for KDE")
            except Exception:
                ax_k.hist(valid_raw, density=True, alpha=0.4, bins=20, label="Density", color=c_raw)
            
            ax_k.axvline(avg, color=c_ma, linestyle="--", linewidth=1.5, label="Mean")
            
            ax_k.set_title(f"{title_prefix} Density Distribution", fontsize=10, fontweight="bold")
            ax_k.set_ylabel("Density", fontweight="bold")
            ax_k.set_xlabel(ylabel, fontweight="bold")
            ax_k.grid(True, alpha=0.3)

            # Legending
            if r == 0 and c == 0:
                ax_ts.set_ylim(top=ax_ts.get_ylim()[1] * 1.5)
                ax_k.set_ylim(top=ax_k.get_ylim()[1] * 1.4)
                for ax_leg in [ax_ts, ax_k]:
                    leg = ax_leg.legend(fontsize=8, framealpha=1.0, loc="upper left", edgecolor="black", fancybox=False)
                    if leg: leg.get_frame().set_linewidth(1.0); leg.set_zorder(10)
    
    fig.text(0.5, 0.02, banner_text, ha="center", va="center", fontsize=10, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none"))

    outroot = OUTROOT_HB if metric == "hbond" else (OUTROOT_CO if metric == "com" else OUTROOT_CT)
    outprefix = os.path.join(outroot, "master_grids", f"MASTER_{title_prefix}_GRID_PUBLICATION_{base.upper()}")
    fig.savefig(f"{outprefix}.png", dpi=300, facecolor="white", format="png", bbox_inches="tight")
    plt.close(fig)
    print(f"[SUCCESS] Saved formatted grid -> {outprefix}.png")

def main():
    ensure_outdirs()
    systems = list_bound_systems()
    if not systems: return print("[WARN] No bound systems found.")

    results = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for fut in as_completed({ex.submit(analyze_one, s): s for s in systems}):
            results.append(fut.result())

    df = build_summary(results)
    if df.empty: return print("[WARN] No data summarized.")
    
    for base, sub_df in df.groupby("base"):
        plot_publication_grid(sub_df, base, metric="hbond")
        plot_publication_grid(sub_df, base, metric="com")
        plot_publication_grid(sub_df, base, metric="contacts")
        
    print("\n[DONE] All metrics analyzed and Publication Master Grids generated successfully.")

if __name__ == "__main__":
    main()
