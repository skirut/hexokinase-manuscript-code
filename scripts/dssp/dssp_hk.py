#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DSSP analysis pipeline for Hexokinase trajectories.

Purpose
-------
Discovers apo and G6P-bound Hexokinase systems, assigns secondary
structure over a configurable trajectory window, caches per-system
results, aggregates secondary-structure statistics, and generates
publication-oriented DSSP plots. The implementation supports an
MDAnalysis DSSP backend and an MDTraj + external DSSP backend.

Inputs
------
- Apo trajectories under ``$HEXOKINASE_DATA_ROOT/structures``
- Bound trajectories under ``$HEXOKINASE_DATA_ROOT/G6P_structure``
- A compatible topology (for example ``structure.prmtop`` or
  ``structure.pdb``)
- A compatible trajectory (for example ``output.xtc``)
- Default analysis window starts at 20 ns and uses pH 3.0, 5.0, and 7.0.

Outputs
-------
Under ``results/dssp/``:
- per-system DSSP assignment matrices and time-series tables
- ``meta.json`` and compressed DSSP caches
- preview plots
- aggregated/statistical tables
- publication plots
- ``logs/dssp_batch.log``

Run from
--------
Run from the repository root with ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/dssp/dssp_hk.py``

Dependencies
------------
Python, NumPy, pandas, MDAnalysis with ``MDAnalysis.analysis.dssp`` or
MDTraj plus an external ``mkdssp``/``dssp`` executable, Matplotlib, and
the repository path utilities.
"""


from __future__ import annotations
import os, sys, json, shutil, time, string
from pathlib import Path
from typing import Optional, Tuple, List, Dict
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings

# Suppress MDAnalysis warnings for cleaner terminal output
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ================= USER SETTINGS =================
from scripts.utils.paths import APO_DIR, BOUND_DIR, DSSP_DIR

ANALYSIS_ROOT = DSSP_DIR

# sampling window & stride
START_NS: float = 20.0           # desired start; will auto-relax if traj shorter
MAX_NS: Optional[float] = None   # cap; None = full after start
STRIDE: Optional[int] = 1        # frame thinning; if None auto-select
TARGET_SAMPLE_FRAMES = 300       # used if STRIDE is None

FALLBACK_DT_PS: float = 10.0     # if times missing

SELECTION: str = "protein"
GUESS_HS: bool = True
BACKEND_PREFERENCE: str = "auto"  # 'auto' | 'mdanalysis' | 'mdtraj'

# MODIFICATION 1: Cap max workers at 4 to prevent Out-Of-Memory crashes
N_WORKERS: Optional[int] = min(4, os.cpu_count() or 4)

MAKE_PREVIEW_PLOTS: bool = True
PUBLISH_PLOTS: bool = True
FORCE_RECOMPUTE: bool = False     # set True to force DSSP recompute even if cache exists

REFERENCE_PH: str = "7.0"
PH_LEVELS = ['3.0', '5.0', '7.0']

# ================= SYSTEM MAPPING =================
SYSTEM_MAPPING = {
    '1qha':       {'apo': '1qha',             'bound': '1QHA_model_with_G6P', 'display': 'HK1'},
    '1qha_G87E':  {'apo': '1qha_G87E',        'bound': '1QHA_G87E_G6P',       'display': 'HK1_G87E'},
    '2nzt':       {'apo': '2nzt',             'bound': '2nzt_model_with_G6P', 'display': 'HK2'},
    '2nzt_D84Y':  {'apo': '2nzt_D84Y',        'bound': '2nzt_D84Y_G6P',       'display': 'HK2_D84Y'},
    '2nzt_D861Y': {'apo': '2nzt_D861Y',       'bound': '2nzt_D861Y_G6P',      'display': 'HK2_D861Y'},
    '2nzt_T680M': {'apo': '2nzt_T680M',       'bound': '2nzt_T680M_G6P',      'display': 'HK2_T680M'},
    '3hm8':       {'apo': '3hm8',             'bound': '3hm8_model_with_G6P', 'display': 'HK3'},
    '3hm8_D867N': {'apo': '3hm8_D867N',       'bound': '3hm8_D867N_G6P',      'display': 'HK3_D867N'},
    '3hm8_D867Y': {'apo': '3hm8_D867Y',       'bound': '3hm8_D867Y_G6P',      'display': 'HK3_D867Y'},
    '3hm8_T546M': {'apo': '3hm8_T546M',       'bound': '3hm8_T546M_G6P',      'display': 'HK3_T546M'},
    '1v4s':       {'apo': '1v4s',             'bound': None,                  'display': 'HK4'},
    '1v4s_D205N': {'apo': '1v4s_D205N',       'bound': None,                  'display': 'HK4_D205N'}
}

# Directory scaffold
PER_SYS_DIR   = ANALYSIS_ROOT / "per_system"
AGG_DIR       = ANALYSIS_ROOT / "aggregated"
STATS_DIR     = ANALYSIS_ROOT / "stats"
PLOTS_DIR     = ANALYSIS_ROOT / "plots"
LOGS_DIR      = ANALYSIS_ROOT / "logs"
for d in [PER_SYS_DIR, AGG_DIR, STATS_DIR, PLOTS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)
LOGFILE = LOGS_DIR / "dssp_batch.log"

def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg.rstrip()}"
    print(line, flush=True)
    try:
        with LOGFILE.open("a") as f:
            f.write(line + "\n")
    except Exception:
        pass

# ---------- discovery ----------
TOPO_PRIORITIES = [
    "structure.prmtop","topol.prmtop","topol.tpr","structure.pdb","protein.pdb",
    "*.prmtop","*.tpr","*.psf","*.pdb",
]
TRAJ_PRIORITIES = ["output.xtc","traj.xtc","*.xtc","*.trr","*.dcd","*.nc","*.crd"]

def _best_match(d: Path, pats: List[str]) -> Optional[Path]:
    for pat in pats:
        hits = sorted(d.glob(pat))
        if hits: return hits[0]
    return None

def find_hk_systems() -> List[Dict]:
    tasks = []
    for sys_key, info in SYSTEM_MAPPING.items():
        display_name = info['display']
        if '_' in display_name:
            enzyme, mutant = display_name.split('_', 1)
        else:
            enzyme, mutant = display_name, 'WT'

        for state in ['apo', 'bound']:
            prefix = info[state]
            if not prefix:
                continue
            base_dir = APO_DIR if state == 'apo' else BOUND_DIR

            for ph in PH_LEVELS:
                sysdir = base_dir / f"{prefix}_pH{ph}" / "prod"
                if not sysdir.exists():
                    sysdir = base_dir / f"{prefix}_pH{ph}"
                if not sysdir.exists():
                    continue

                top = _best_match(sysdir, TOPO_PRIORITIES)
                traj = _best_match(sysdir, TRAJ_PRIORITIES)
                if top and traj:
                    state_display = state.capitalize()
                    system_id = f"{display_name}_{state_display}_pH{ph}"
                    
                    tasks.append({
                        'sysdir': sysdir,
                        'top': top,
                        'traj': traj,
                        'system_id': system_id,
                        'isoform': f"{display_name}_{state_display}",
                        'enzyme': enzyme,
                        'mutant': mutant,
                        'state': state_display,
                        'pH': ph,
                        'replicate': 'rep1'
                    })
    return tasks

# ---------- backends ----------
def find_dssp_exec() -> Optional[str]:
    return shutil.which("mkdssp") or shutil.which("dssp")

def choose_backend(pref: str = "auto") -> str:
    def have_mda():
        try:
            import MDAnalysis  # noqa
            from MDAnalysis.analysis import dssp as _dssp  # noqa
            return True
        except Exception:
            return False
    def have_mdtraj():
        try:
            import mdtraj  # noqa
            return find_dssp_exec() is not None
        except Exception:
            return False
    if pref == "mdanalysis":
        if have_mda():
            log("[BACKEND] Using MDAnalysis.analysis.dssp")
            return "mdanalysis"
        raise RuntimeError("Requested mdanalysis backend unavailable.")
    if pref == "mdtraj":
        if have_mdtraj():
            log(f"[BACKEND] Using MDTraj + {find_dssp_exec()}")
            return "mdtraj"
        raise RuntimeError("Requested mdtraj backend unavailable.")
    if have_mdtraj():
        log(f"[BACKEND] Auto-select: MDTraj + {find_dssp_exec()}")
        return "mdtraj"
    if have_mda():
        log("[BACKEND] Auto-select: MDAnalysis")
        return "mdanalysis"
    raise RuntimeError("No DSSP backend available (install mdtraj+mkdssp or MDAnalysis).")

# ---------- helpers ----------
def _auto_relax_start_ns(times_ns: np.ndarray, want_start: float) -> float:
    if len(times_ns)==0: return 0.0
    tmax = float(np.nanmax(times_ns))
    if tmax < want_start:
        return 0.0
    return want_start

def _auto_stride_for_target(n_frames: int, target: int = TARGET_SAMPLE_FRAMES) -> int:
    if n_frames <= target or target <= 0:
        return 1
    return max(1, int(np.ceil(n_frames / target)))

# ---------- MDTraj runner (CHUNKED FOR MEMORY SAFETY) ----------
def run_dssp_mdtraj(top: Path, traj: Path) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    import mdtraj as md, os as _os
    exe = find_dssp_exec()
    if exe: _os.environ["DSSP"] = exe

    # Fallback to alternative topologies if PDB/PRMTOP fail
    top_to_use = top
    if top.name.endswith(".prmtop"):
        cand = top.with_name("structure.pdb")
        if cand.exists(): top_to_use = cand
    elif top.name.endswith(".pdb"):
        cand = top.with_name("structure.prmtop")
        if cand.exists(): top_to_use = cand

    CHUNK_SIZE = 500
    times_all_list = []
    ss8_list = []
    resids = None

    try:
        # Load and compute in chunks to avoid OOM killer
        for chunk in md.iterload(str(traj), top=str(top_to_use), chunk=CHUNK_SIZE):
            sel = chunk.topology.select("protein")
            if sel.size == 0:
                raise RuntimeError("MDTraj selection 'protein' empty.")
            chunk = chunk.atom_slice(sel)
            
            if resids is None:
                resids = np.array([r.resSeq for r in chunk.topology.residues])
            
            if hasattr(chunk, "time") and chunk.time is not None and len(chunk.time) > 0:
                times_all_list.extend(chunk.time / 1000.0)
            else:
                # Estimate missing times
                start_idx = len(times_all_list)
                times_all_list.extend(np.arange(start_idx, start_idx + chunk.n_frames) * (FALLBACK_DT_PS/1000.0))

            chunk_ss8 = md.compute_dssp(chunk, simplified=False)
            ss8_list.append(chunk_ss8)
    except Exception as e:
        raise RuntimeError(f"MDTraj chunked processing failed for {top_to_use} + {traj}: {e}")

    if not ss8_list:
        raise RuntimeError(f"No frames loaded or processed for {traj}.")

    times_all = np.array(times_all_list)
    ss8 = np.concatenate(ss8_list, axis=0)

    use_start = _auto_relax_start_ns(times_all, START_NS)
    mask = times_all >= use_start
    if MAX_NS is not None:
        mask &= (times_all <= MAX_NS)
    if not mask.any():
        raise RuntimeError("No frames found in selected window (MDTraj).")

    if STRIDE and STRIDE > 1:
        idx = np.where(mask)[0][::STRIDE]
    else:
        idx_all = np.where(mask)[0]
        stride_here = _auto_stride_for_target(len(idx_all), TARGET_SAMPLE_FRAMES) if STRIDE is None else 1
        idx = idx_all[::stride_here]

    times_ns = times_all[idx]
    ss8_filtered = ss8[idx]

    map3 = {'H':'H','G':'H','I':'H','E':'E','B':'E'}
    labels3 = np.array([[map3.get(c, 'C') for c in row] for row in ss8_filtered])
    raw = np.array([['-' if c not in ('H','E','G','I','B') else ('H' if c in ('H','G','I') else 'E')
                     for c in row] for row in ss8_filtered])
    return times_ns, resids, labels3, raw

# ---------- MDAnalysis runner ----------
def run_dssp_mdanalysis(top: Path, traj: Path) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP

    def _try_load(topfile: Path, trajfile: Path):
        return mda.Universe(str(topfile), str(trajfile))

    u = None
    errA = None
    try:
        u = _try_load(top, traj)
    except Exception as e:
        errA = e
        alt_top = None
        if top.name.endswith(".prmtop"):
            cand = top.with_name("structure.pdb")
            if cand.exists(): alt_top = cand
        elif top.name.endswith(".pdb"):
            cand = top.with_name("structure.prmtop")
            if cand.exists(): alt_top = cand
        if alt_top:
            try:
                u = _try_load(alt_top, traj)
                log(f"[INFO] Loaded with alternate topology: {alt_top}")
            except Exception as e2:
                errA = e2
    if u is None:
        raise RuntimeError(f"MDAnalysis load failed for {top} + {traj}: {errA}")

    ag = u.select_atoms(SELECTION)
    if ag.n_atoms == 0:
        raise RuntimeError(f"Selection yielded 0 atoms: {SELECTION}")

    # HK Bulletproof filter for valid backbone atoms
    if SELECTION == 'protein':
        valid_residues = []
        for res in ag.residues:
            names = list(res.atoms.names)
            if (names.count('N') == 1 and names.count('CA') == 1 and 
                names.count('C') == 1 and names.count('O') == 1):
                valid_residues.append(res)
        ag = mda.core.groups.ResidueGroup(valid_residues).atoms

    times_all = []
    frames_all = []
    for ts in u.trajectory:
        t_ps = ts.time if hasattr(ts, "time") else None
        t_ns = (t_ps/1000.0) if (t_ps is not None and not np.isnan(t_ps)) else None
        if t_ns is None:
            t_ns = ts.frame * (FALLBACK_DT_PS/1000.0)
        times_all.append(t_ns)
        frames_all.append(int(ts.frame))
    times_all = np.asarray(times_all, dtype=float)
    frames_all = np.asarray(frames_all, dtype=int)

    use_start = _auto_relax_start_ns(times_all, START_NS)
    mask = (times_all >= use_start)
    if MAX_NS is not None:
        mask &= (times_all <= MAX_NS)
    if not mask.any():
        raise RuntimeError("No frames found in selected window (MDAnalysis).")

    if STRIDE and STRIDE>1:
        idxs = np.where(mask)[0][::STRIDE]
    else:
        idx_all = np.where(mask)[0]
        stride_here = _auto_stride_for_target(len(idx_all), TARGET_SAMPLE_FRAMES) if STRIDE is None else 1
        idxs = idx_all[::stride_here]

    times_ns = times_all[idxs]
    target_frames = frames_all[idxs].tolist()

    dssp = DSSP(ag, guess_hydrogens=GUESS_HS, update_selections=False)
    dssp.run()
    full_dssp = dssp.results.dssp
    resids = np.array(dssp.results.resids)

    try:
        if full_dssp.shape[0] > max(target_frames):
            labels = full_dssp[target_frames, :]
        elif full_dssp.shape[0] == len(times_all):
            labels = full_dssp[idxs, :]
        else:
            labels = full_dssp[idxs, :]
    except Exception:
        labels = full_dssp[idxs, :]

    labels3 = np.where(labels == '-', 'C', labels)
    raw = labels
    return np.array(times_ns), resids, labels3, raw

# ---------- outputs ----------
def write_per_system(task: Dict, times_ns, resids, labels3, raw, outdir: Optional[Path] = None):
    system_id = task['system_id']
    if outdir is None:
        outdir = PER_SYS_DIR / task['state'] / system_id
    outdir.mkdir(parents=True, exist_ok=True)
    
    df_mat = pd.DataFrame(labels3.T, index=resids,
                          columns=[f"{t:.3f}ns" for t in times_ns])
    df_mat.index.name="resid"
    df_mat.to_csv(outdir/"dssp_assignments_matrix.tsv", sep="\t")

    frac_H = (labels3=='H').mean(axis=1)
    frac_E = (labels3=='E').mean(axis=1)
    frac_C = (labels3=='C').mean(axis=1)
    pd.DataFrame({"time_ns":times_ns,"H":frac_H,"E":frac_E,"C":frac_C}) \
        .to_csv(outdir/"dssp_fraction_timeseries.csv", index=False)

    prop_H = (labels3=='H').mean(axis=0)
    prop_E = (labels3=='E').mean(axis=0)
    prop_C = (labels3=='C').mean(axis=0)
    pd.DataFrame({"resid":resids,"H":prop_H,"E":prop_E,"C":prop_C}) \
        .to_csv(outdir/"dssp_per_residue_propensity.csv", index=False)

    meta_data = task.copy()
    meta_data.update({
        "n_frames": int(labels3.shape[0]), 
        "n_residues": int(labels3.shape[1]),
        "start_ns_requested": START_NS, 
        "start_ns_effective": float(times_ns[0]) if len(times_ns)>0 else 0.0,
        "max_ns": MAX_NS, "stride": STRIDE, "selection": SELECTION,
        "top": str(task['top']), "traj": str(task['traj'])
    })
    meta_data.pop('sysdir', None)
    
    (outdir/"meta.json").write_text(json.dumps(meta_data, indent=2))

    try:
        np.savez_compressed(outdir/"dssp.npz", times_ns=times_ns, resids=resids, labels3=labels3, raw=raw)
    except Exception:
        log(f"[WARN] Failed to write cache for {system_id}")

    if MAKE_PREVIEW_PLOTS:
        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import ListedColormap
            mpl = __import__("matplotlib")
            mpl.rcParams.update({
                "font.family": "sans-serif",
                "font.weight": "bold",
                "axes.labelweight": "bold",
                "axes.titleweight": "bold",
                "font.size": 10,
                "axes.titlesize": 11,
                "axes.labelsize": 10,
                "figure.dpi": 1200,
                "savefig.dpi": 1200,
            })

            Z = np.vectorize({'H':0,'E':1,'C':2}.get)(labels3).T
            cmap = ListedColormap(["#1f77b4","#ff7f0e","#7f7f7f"])
            plt.figure(figsize=(8,4), dpi=1200)
            plt.imshow(Z, aspect='auto', origin='lower', interpolation='nearest', cmap=cmap)
            plt.xlabel("Frame", fontweight='bold'); plt.ylabel("Residue", fontweight='bold')
            plt.xticks(fontweight='bold')
            plt.yticks(fontweight='bold')
            nx = Z.shape[1]
            xt = np.linspace(0, nx-1, min(6, nx), dtype=int)
            plt.xticks(xt, [f"{times_ns[i]:.1f}" for i in xt], fontweight='bold')
            cbar = plt.colorbar(ticks=[0,1,2]); cbar.set_ticklabels(["Helix","Sheet","Coil"])
            plt.tight_layout(); plt.savefig(outdir/"preview_timeline.png"); plt.close()

            F = pd.DataFrame({"H":frac_H,"E":frac_E,"C":frac_C}, index=times_ns)
            plt.figure(figsize=(6, 3.2), dpi=1200)
            for k, col in zip(["H","E","C"], ["#1f77b4","#ff7f0e","#7f7f7f"]):
                plt.plot(F.index.values, F[k].values, label=k, color=col, linewidth=2.5)
            plt.xlabel("Time (ns)", fontsize=11, fontweight='bold')
            plt.ylabel("Fraction", fontsize=11, fontweight='bold')
            plt.xticks(fontweight='bold', fontsize=10)
            plt.yticks(fontweight='bold', fontsize=10)
            plt.tight_layout(); plt.savefig(outdir/"preview_fractions.png"); plt.close()
        except Exception as e:
            log(f"[WARN] preview plotting failed: {e}")

# ---------- worker ----------
def process_system(task: Dict, backend: str):
    system_id = task['system_id']
    top = task['top']
    traj = task['traj']
    
    outdir = PER_SYS_DIR / task['state'] / system_id
    outdir.mkdir(parents=True, exist_ok=True)

    cachef = outdir/"dssp.npz"
    if cachef.exists() and not FORCE_RECOMPUTE:
        try:
            data = np.load(cachef, allow_pickle=True)
            times_ns = data["times_ns"]
            resids = data["resids"]
            labels3 = data["labels3"]
            raw = data["raw"]
            write_per_system(task, times_ns, resids, labels3, raw, outdir)
            return system_id, True, f"OK (cache)"
        except Exception as e:
            log(f"[WARN] Failed to load cache {cachef}: {e}")

    try_order = [(top, traj)]
    if top.name.endswith(".pdb"):
        alt = top.with_name("structure.prmtop")
        if alt.exists():
            try_order.insert(0, (alt, traj))
    elif top.name.endswith(".prmtop"):
        alt = top.with_name("structure.pdb")
        if alt.exists():
            try_order.append((alt, traj))

    last_err = None
    start_t = time.time()
    for ttop, ttraj in try_order:
        try:
            if backend == "mdanalysis":
                times_ns, resids, labels3, raw = run_dssp_mdanalysis(ttop, ttraj)
            else:
                times_ns, resids, labels3, raw = run_dssp_mdtraj(ttop, ttraj)
            write_per_system(task, times_ns, resids, labels3, raw, outdir)
            dur = time.time() - start_t
            return system_id, True, f"OK ({len(times_ns)} frames, {len(resids)} residues, {dur:.1f}s)"
        except Exception as e:
            last_err = e
            log(f"[WARN] DSSP attempt failed for {ttop}:{ttraj} -> {e}")

    if backend == "mdanalysis":
        try:
            times_ns, resids, labels3, raw = run_dssp_mdtraj(top, traj)
            write_per_system(task, times_ns, resids, labels3, raw, outdir)
            return system_id, True, f"OK via MDTraj fallback ({len(times_ns)} frames)"
        except Exception as e2:
            last_err = (last_err, e2)
            log(f"[WARN] fallback MDTraj also failed: {e2}")

    return system_id, False, f"{last_err}"

# ---------- aggregation & stats ----------
def aggregate_per_system():
    rows_global, rows_res = [], []
    
    for state_dir in PER_SYS_DIR.iterdir():
        if not state_dir.is_dir(): continue
        for sysdir in state_dir.iterdir():
            if not sysdir.is_dir(): continue
            meta_f = sysdir/"meta.json"
            frac_f = sysdir/"dssp_fraction_timeseries.csv"
            prop_f = sysdir/"dssp_per_residue_propensity.csv"
            if not (meta_f.exists() and frac_f.exists() and prop_f.exists()):
                continue
            meta = json.loads(meta_f.read_text())
            dfF = pd.read_csv(frac_f)
            rows_global.append({
                "system_id": meta.get("system_id", "Unknown"), "isoform": meta.get("isoform", "Unknown"),
                "pH": str(meta.get("pH", "Unknown")), "replicate": meta.get("replicate", "Unknown"),
                "n_frames": meta.get("n_frames", int(dfF.shape[0])),
                "mean_H": dfF["H"].mean(), "mean_E": dfF["E"].mean(), "mean_C": dfF["C"].mean()
            })
            dfP = pd.read_csv(prop_f)
            dfP.insert(0,"system_id",meta.get("system_id", "Unknown"))
            dfP.insert(1,"isoform",meta.get("isoform", "Unknown"))
            dfP.insert(2,"pH",str(meta.get("pH", "Unknown")))
            dfP.insert(3,"replicate",meta.get("replicate", "Unknown"))
            rows_res.append(dfP)
            
    G = pd.DataFrame(rows_global)
    R = pd.concat(rows_res, ignore_index=True) if rows_res else pd.DataFrame(
        columns=["system_id","isoform","pH","replicate","resid","H","E","C"]
    )
    AGG_DIR.mkdir(parents=True, exist_ok=True)
    G.to_csv(AGG_DIR/"global_mean_fraction_per_system.csv", index=False)
    R.to_csv(AGG_DIR/"per_residue_propensity_per_system.csv", index=False)
    return G, R

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return np.minimum(out, 1.0)

def kw_mwu_within_isoform(G: pd.DataFrame, metric: str, prefix: str):
    from scipy.stats import kruskal, mannwhitneyu
    kw_rows, mwu_rows = [], []
    for iso, dfI in G.groupby("isoform"):
        groups = [(ph, sub[metric].dropna().values) for ph, sub in dfI.groupby("pH")]
        groups = [(ph, v) for ph, v in groups if v.size>0]
        if len(groups) < 2: continue
        try: stat, p = kruskal(*[v for _,v in groups], nan_policy='omit')
        except Exception: stat, p = np.nan, np.nan
        kw_rows.append({"isoform":iso, "metric":metric, "kw_stat":stat, "kw_p":p, "n_groups":len(groups)})
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                g1,v1 = groups[i]; g2,v2 = groups[j]
                if v1.size==0 or v2.size==0: continue
                try: u,p2 = mannwhitneyu(v1, v2, alternative='two-sided')
                except Exception: u,p2 = np.nan, np.nan
                mwu_rows.append({"isoform":iso,"metric":metric,"pH1":g1,"pH2":g2,"u_stat":u,"p_raw":p2,
                                 "n1":int(v1.size),"n2":int(v2.size)})
    df_kw = pd.DataFrame(kw_rows)
    if not df_kw.empty:
        df_kw["kw_q_bh"] = bh_fdr(df_kw["kw_p"].fillna(1.0).values)
        df_kw.to_csv(STATS_DIR/f"{prefix}_kw_omnibus.csv", index=False)
    df_mwu = pd.DataFrame(mwu_rows)
    if not df_mwu.empty:
        df_mwu["q_bh"] = bh_fdr(df_mwu["p_raw"].fillna(1.0).values)
        df_mwu.to_csv(STATS_DIR/f"{prefix}_mwu_pairwise.csv", index=False)

def anova_tukey_within_isoform(G: pd.DataFrame, metric: str, prefix: str):
    try:
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
    except Exception:
        log("[INFO] statsmodels not found; skipping ANOVA/Tukey.")
        return
    rows_aov, rows_tk = [], []
    for iso, dfI in G.groupby("isoform"):
        if dfI["pH"].nunique() < 2: continue
        try:
            m = ols(f"{metric} ~ C(pH)", data=dfI).fit()
            aov = sm.stats.anova_lm(m, typ=2)
            p = float(aov["PR(>F)"].iloc[0])
            rows_aov.append({"isoform":iso,"metric":metric,"anova_p":p})
            tk = pairwise_tukeyhsd(dfI[metric], dfI["pH"], alpha=0.05)
            df_tk = pd.DataFrame(tk._results_table.data[1:], columns=tk._results_table.data[0])
            df_tk.insert(0,"isoform",iso); df_tk.insert(1,"metric",metric)
            rows_tk.append(df_tk)
        except Exception:
            continue
    if rows_aov:
        df_aov = pd.DataFrame(rows_aov)
        df_aov["anova_q_bh"] = bh_fdr(df_aov["anova_p"].fillna(1.0).values)
        df_aov.to_csv(STATS_DIR/f"{prefix}_anova.csv", index=False)
    if rows_tk:
        pd.concat(rows_tk, ignore_index=True).to_csv(STATS_DIR/f"{prefix}_tukey.csv", index=False)

def residue_level_stats(R: pd.DataFrame, state: str):
    from scipy.stats import kruskal, mannwhitneyu
    kw_rows, mwu_rows = [], []
    for (iso, resid), dfIR in R.groupby(["isoform","resid"]):
        groups = [(ph, sub[state].dropna().values) for ph, sub in dfIR.groupby("pH")]
        groups = [(ph,v) for ph,v in groups if v.size>0]
        if len(groups)<2: continue
        try: stat,p = kruskal(*[v for _,v in groups], nan_policy='omit')
        except Exception: stat,p = np.nan, np.nan
        kw_rows.append({"isoform":iso,"resid":resid,"state":state,"kw_stat":stat,"kw_p":p,"n_groups":len(groups)})
        for i in range(len(groups)):
            for j in range(i+1,len(groups)):
                g1,v1 = groups[i]; g2,v2 = groups[j]
                try: u,p2 = mannwhitneyu(v1,v2,alternative='two-sided')
                except Exception: u,p2 = np.nan, np.nan
                mwu_rows.append({"isoform":iso,"resid":resid,"state":state,"pH1":g1,"pH2":g2,"u_stat":u,"p_raw":p2,
                                 "n1":int(v1.size),"n2":int(v2.size)})
    if kw_rows:
        df = pd.DataFrame(kw_rows); df["kw_q_bh"] = bh_fdr(df["kw_p"].fillna(1.0).values)
        df.to_csv(STATS_DIR/f"residue_{state}_kw.csv", index=False)
    if mwu_rows:
        df = pd.DataFrame(mwu_rows); df["q_bh"] = bh_fdr(df["p_raw"].fillna(1.0).values)
        df.to_csv(STATS_DIR/f"residue_{state}_mwu.csv", index=False)

# ---------- publication composite timeline grid generation ----------
def generate_publication_composite_grids():
    log("[PLOTS] Generating publication composite timeline grids (Direct Plotting)...")
    
    # Gather data from computed system directories
    records = []
    for state_dir in PER_SYS_DIR.iterdir():
        if not state_dir.is_dir(): continue
        for sdir in state_dir.iterdir():
            if not sdir.is_dir(): continue
            
            meta_f = sdir / "meta.json"
            csv_f = sdir / "dssp_fraction_timeseries.csv"
            if not meta_f.exists() or not csv_f.exists(): continue

            try:
                meta = json.loads(meta_f.read_text())
                records.append({
                    "enzyme": meta.get("enzyme", "HK"),
                    "mut": meta.get("mutant", "WT"),
                    "state": meta.get("state", "Apo"),
                    "ph": str(meta.get("pH", "")),
                    "csv_f": csv_f
                })
            except Exception:
                continue

    if not records:
        log("[WARN] No valid system records found for composite grids.")
        return

    grouped = {}
    for r in records:
        grouped.setdefault(r["enzyme"], []).append(r)

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import pandas as pd

    def sort_conditions(cond):
        mut, state = cond
        mut_score = 0 if mut == 'WT' else 1
        state_score = 0 if state.lower() == 'apo' else 1
        return (mut_score, mut, state_score)

    for enz, recs in grouped.items():
        conditions = list(set((r['mut'], r['state']) for r in recs))
        sorted_conditions = sorted(conditions, key=sort_conditions)

        nrows = len(sorted_conditions)
        ncols = len(PH_LEVELS)
        if nrows == 0:
            continue

        # --- DYNAMIC Y-AXIS CALCULATION ---
        system_max_frac = 0.0
        for r in recs:
            if r['csv_f'].exists():
                df_temp = pd.read_csv(r['csv_f'])
                max_in_file = df_temp[['H', 'E', 'C']].max().max()
                if max_in_file > system_max_frac:
                    system_max_frac = max_in_file
                    
        y_lim_upper = min(1.0, np.ceil((system_max_frac + 0.08) * 20) / 20.0)

        fig_w = 6.0 * ncols
        fig_h = 3.5 * nrows + 1.5

        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharex=True, sharey=True)

        if nrows == 1 and ncols == 1:
            axes = [[axes]]
        elif nrows == 1:
            axes = [axes]
        elif ncols == 1:
            axes = [[ax] for ax in axes]
        
        fig.suptitle(f"DSSP Analysis for {enz.upper()}", fontsize=26, fontweight="bold", y=0.98)
        letters = string.ascii_uppercase

        for row_idx, cond in enumerate(sorted_conditions):
            mut, state = cond
            row_letter = letters[row_idx] if row_idx < len(letters) else f"R{row_idx+1}"

            for col_idx, ph_val in enumerate(PH_LEVELS):
                ax = axes[row_idx][col_idx]

                matched_file = None
                for r in recs:
                    if r['mut'] == mut and r['state'] == state and r['ph'] == ph_val:
                        matched_file = r['csv_f']
                        break

                if matched_file and matched_file.exists():
                    df = pd.read_csv(matched_file)
                    
                    ax.plot(df['time_ns'], df['H'], color='#1f77b4', lw=2.5, label='Alpha-Helix')
                    ax.plot(df['time_ns'], df['E'], color='#ff7f0e', lw=2.5, label='Beta-Sheet')
                    ax.plot(df['time_ns'], df['C'], color='#7f7f7f', lw=2.5, label='Coil / Turn / Other')
                    
                    ax.tick_params(axis='both', which='major', labelsize=14)
                    for spine in ax.spines.values():
                        spine.set_linewidth(1.5)
                        spine.set_color('#333333')
                        
                    ax.set_ylim(0.0, y_lim_upper)
                else:
                    ax.text(0.5, 0.5, "Data Missing", ha='center', va='center', fontsize=16, color='gray')

                if row_idx == 0:
                    ax.set_title(f"pH {ph_val}", fontsize=20, fontweight="bold", pad=15)

                if col_idx == ncols - 1:
                    row_label = f"{mut} - {state}"
                    ax.text(1.05, 0.5, row_label, transform=ax.transAxes, 
                            fontsize=20, fontweight="bold", ha="left", va="center", rotation=-90)

                if col_idx == 0:
                    ax.text(-0.15, 0.5, f"{row_letter})", transform=ax.transAxes, 
                            fontsize=20, fontweight="bold", ha="right", va="center")

        fig.subplots_adjust(left=0.08, right=0.92, top=0.88, bottom=0.15, wspace=0.1, hspace=0.1)

        fig.supxlabel("Time (ns)", fontsize=22, fontweight='bold', y=0.08)
        fig.supylabel("Fraction", fontsize=22, fontweight='bold', x=0.02)

        legend_elements = [
            Line2D([0], [0], color='#1f77b4', lw=4.0, label='Alpha-Helix Fraction'),
            Line2D([0], [0], color='#ff7f0e', lw=4.0, label='Beta-Sheet Fraction'),
            Line2D([0], [0], color='#7f7f7f', lw=4.0, label='Coil / Turn / Other'),
        ]

        fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.01),
                   ncol=3, frameon=True, facecolor='white', edgecolor='black', 
                   prop={'weight': 'bold', 'size': 18})

        out_png = f"{enz}_composite_preview_fractions_PUBLICATION.png"
        out_pdf = f"{enz}_composite_preview_fractions_PUBLICATION.pdf"

        fig.savefig(PLOTS_DIR / out_png, dpi=1200, bbox_inches='tight', facecolor='white')
        fig.savefig(PLOTS_DIR / out_pdf, dpi=1200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        log(f"[PLOTS] Composite generated -> {out_png}")

# ---------- standard plots ----------
def publish_plots(G: pd.DataFrame, R: pd.DataFrame):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.dpi": 1200,
        "savefig.dpi": 1200,
    })

    for iso, dfI in G.groupby("isoform"):
        phs = sorted(dfI["pH"].unique(), key=lambda x: float(x) if str(x).replace('.','',1).isdigit() else 1e9)
        x = np.arange(len(phs))
        means, sems = [], []
        for ph in phs:
            vals = dfI.loc[dfI["pH"]==ph, "mean_H"].dropna().values
            means.append(vals.mean() if vals.size else np.nan)
            sems.append(vals.std(ddof=1)/np.sqrt(vals.size) if vals.size>1 else 0.0)
        plt.figure(dpi=1200, figsize=(6,4))
        plt.bar(x, means, yerr=sems, capsize=3)
        for i, ph in enumerate(phs):
            pts = dfI.loc[dfI["pH"]==ph, "mean_H"].dropna().values
            if pts.size:
                jitter = (np.random.rand(pts.size)-0.5)*0.15
                plt.scatter(np.full(pts.shape,i)+jitter, pts, s=20)
        plt.xticks(x, phs, fontweight='bold')
        plt.yticks(fontweight='bold')
        plt.ylabel("Mean helix fraction", fontweight='bold')
        plt.title(f"{iso}: helicity vs pH", fontweight='bold')
        plt.tight_layout(); plt.savefig(PLOTS_DIR/f"{iso}_global_helix_vs_pH.png"); plt.close()

    ref = REFERENCE_PH
    for iso, dfI in R.groupby("isoform"):
        pivot = dfI.pivot_table(index="resid", columns="pH", values="H", aggfunc="mean")
        # Ensure we operate on percent vs fraction based on data, converting to % for delta plot
        pivot = pivot * 100.0 
        cols_map = {str(c): c for c in pivot.columns}
        if ref not in cols_map:
            continue
        ref_col = cols_map[ref]
        delta = pivot.subtract(pivot[ref_col], axis=0)
        
        plt.figure(dpi=1200, figsize=(10,4))
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)
        im = plt.imshow(delta.values, aspect='auto', origin='lower', interpolation='nearest', cmap='RdBu_r', norm=norm)
        cbar = plt.colorbar(im)
        cbar.set_label(f"Δ Helicity % (vs pH {ref})", fontweight='bold')
        plt.xlabel("pH", fontweight='bold')
        plt.ylabel("Residue", fontweight='bold')
        plt.xticks(np.arange(len(delta.columns)), [str(c) for c in delta.columns], rotation=45, fontweight='bold')
        nres = delta.shape[0]
        yt = np.linspace(0, nres-1, min(10, nres), dtype=int)
        plt.yticks(yt, delta.index[yt], fontweight='bold')
        plt.title(f"{iso}: Δ-Helicity per residue", fontweight='bold')
        plt.tight_layout(); plt.savefig(PLOTS_DIR/f"{iso}_delta_helix_vs_pH{ref}.png"); plt.close()

    # Generate publication composite grids directly
    generate_publication_composite_grids()

# ---------- main ----------
def main():
    log("# ===== HK DSSP PUBLICATION PIPELINE (parallel) =====")
    log(f"Output root: {ANALYSIS_ROOT}")

    tasks = find_hk_systems()
    if not tasks:
        log("[ERROR] No systems found."); sys.exit(1)

    try:
        backend = choose_backend(BACKEND_PREFERENCE)
    except Exception as e:
        log(str(e)); sys.exit(2)

    log(f"Discovered {len(tasks)} tasks. Using {N_WORKERS} workers.")

    ok = 0
    with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(tasks))) as ex:
        futs = {ex.submit(process_system, task, backend): task for task in tasks}
        for fut in as_completed(futs):
            task = futs[fut]
            try:
                system_id, success, msg = fut.result()
                if success:
                    ok += 1
                    log(f"[OK] {system_id}: {msg}")
                else:
                    log(f"[FAIL] {system_id}: {msg}")
            except Exception as e:
                log(f"[FAIL] {task['system_id']}: {e}")

    log(f"[DONE] Per-system stage: {ok}/{len(tasks)} succeeded.")

    G, R = aggregate_per_system()
    log(f"[AGG] Global rows: {len(G)} ; Residue rows: {len(R)}")

    for metric in ["mean_H","mean_E","mean_C"]:
        kw_mwu_within_isoform(G, metric, prefix=f"global_{metric}")
        anova_tukey_within_isoform(G, metric, prefix=f"global_{metric}")

    if not R.empty:
        for state in ["H","E","C"]:
            residue_level_stats(R, state)

    if PUBLISH_PLOTS:
        publish_plots(G, R)

    log("# ===== ALL COMPLETED SUCCESSFULLY =====")
    log(f"Per-system -> {PER_SYS_DIR}\nAggregated -> {AGG_DIR}\nStats -> {STATS_DIR}\nPlots -> {PLOTS_DIR}")

if __name__ == "__main__":
    main()
