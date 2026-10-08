#!/usr/bin/env python3
"""Comprehensive Hexokinase surface, propensity, ASA, and variable-site analysis.

Purpose
-------
Parses full Hexokinase PDB files, MolLayers surface-layer PDB files, and
NACCESS RSA files; calculates surface composition, log propensity, ASA,
and variable-site statistics; and generates the publication figures
implemented in the script.

Inputs
------
- MolLayers surface/layer files supplied through ``--layers``.
- NACCESS ``.rsa`` files supplied through ``--rsa``.
- Output directory supplied through ``--outdir``.
- The input filenames and residue labels must match the parsing logic
  in the script.

Outputs
-------
Four publication figures in the requested output directory:
- ``hk_surface_composition_4panel.png``
- ``hk_propensity_asa_2panel.png``
- ``hk_variable_sites_heatmap.png``
- ``hk_comprehensive_mega_plot.png``

Run from
--------
Run from the repository root or another directory where the script is
importable; this script uses explicit command-line input directories.

Usage
-----
``python scripts/surface_analysis/hk_layers_analysis.py --layers <layers_dir> --rsa <rsa_dir> --outdir <output_dir>``

Dependencies
------------
Python, NumPy, pandas, SciPy, Matplotlib, and the supplied MolLayers
and NACCESS-derived input files.
"""


import os
import csv
import argparse
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory as btf
from scipy.stats import chisquare, fisher_exact

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------
ISO_LABELS = {
    "1qha": "HK1 (1QHA)",
    "2nzt": "HK2 (2NZT)",
    "3hm8": "HK3 (3HM8)",
    "1v4s": "HK4 (1V4S)"
}

TITRATABLE = {"ASP", "GLU", "HIS", "HID", "HIE", "HIP", "LYS", "ARG"}
PSEUDOCOUNT = 1e-6

AA_CLASSES = {
    'ARG': 'Positively Charged', 'LYS': 'Positively Charged', 'HIS': 'Positively Charged', 
    'HID': 'Positively Charged', 'HIE': 'Positively Charged', 'HIP': 'Positively Charged',
    'ASP': 'Negatively Charged', 'GLU': 'Negatively Charged',
    'SER': 'Polar', 'THR': 'Polar', 'ASN': 'Polar', 'GLN': 'Polar', 'TYR': 'Polar', 'CYS': 'Polar',
    'GLY': 'Nonpolar', 'ALA': 'Nonpolar', 'VAL': 'Nonpolar', 'LEU': 'Nonpolar', 'ILE': 'Nonpolar', 
    'MET': 'Nonpolar', 'PRO': 'Nonpolar', 'PHE': 'Nonpolar', 'TRP': 'Nonpolar'
}
CLASS_LABELS = ["Nonpolar", "Polar", "Positively Charged", "Negatively Charged"]

AA_3_TO_1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'HID': 'H', 'HIE': 'H', 'HIP': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T',
    'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Color-Blind Friendly Palette (Okabe-Ito inspired)
COLOR_PROTEIN       = "#0072B2"  # Blue
COLOR_SURFACE       = "#D55E00"  # Vermillion
COLOR_PROTEIN_LIGHT = "#56B4E9"  # Sky Blue
COLOR_SURFACE_LIGHT = "#E69F00"  # Orange
COLOR_MAP_PROP      = "coolwarm" # Diverging colorblind friendly

# Composition Colors
CLASS_COLORS = {
    "Nonpolar": "#009E73",           # Bluish Green
    "Polar": "#F0E442",              # Yellow
    "Positively Charged": "#0072B2", # Blue
    "Negatively Charged": "#D55E00", # Vermillion
}

# Heatmap Swatch Colors
HEATMAP_COLORS = {
    "basic":  "#0072B2",   # Blue
    "acidic": "#D55E00",   # Vermillion
    "other":  "#CC79A7",   # Reddish Purple
}

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 15,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# ------------------------------------------------------------
# CORE PARSING & MATH FUNCTIONS
# ------------------------------------------------------------
def tit_group_letter(letter: str):
    if letter in ("D", "E"): return "acidic"
    if letter in ("H", "K", "R"): return "basic"
    if letter in ("Y", "C"): return "other"
    return None

def get_csv_val(r, iso, suffix):
    if f"{iso}_{suffix}" in r: return r[f"{iso}_{suffix}"]
    return ""

def parse_pdb_stats(pdb_path):
    atoms, residues, resnum_only, resnum_to_name = 0, set(), set(), {}
    if not os.path.exists(pdb_path): return atoms, residues, resnum_only, resnum_to_name
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM "):
                atoms += 1
                resname, chain, resnum = line[17:20].strip().upper(), line[21].strip(), line[22:26].strip()
                residues.add((resname, chain, resnum))
                resnum_only.add(resnum)
                resnum_to_name[resnum] = resname
    return atoms, residues, resnum_only, resnum_to_name

def parse_rsa(path):
    out = {}
    if not os.path.exists(path): return out
    with open(path) as fh:
        for raw in fh:
            if not raw.strip() or not raw.lstrip().startswith("RES"): continue
            line = raw.rstrip("\n")
            if len(line) < 20: continue
            try:
                resname, resnum = line[4:7].strip().upper(), line[9:13].strip()
                abs_val = float(line[13:22].strip())
                out[resnum] = (resname, abs_val)
            except Exception:
                continue
    return out

def get_composition(res_set):
    counts = {c: 0 for c in CLASS_LABELS}
    valid_res = 0
    for (resname, _, _) in res_set:
        c = AA_CLASSES.get(resname, None)
        if c is not None:
            counts[c] += 1
            valid_res += 1
    if valid_res == 0: return [0]*4, [0]*4
    pcts = [(counts[c] / valid_res) * 100 for c in CLASS_LABELS]
    abs_counts = [counts[c] for c in CLASS_LABELS]
    return pcts, abs_counts

def bh_fdr(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    if n == 0: return pvals
    sort_idx = np.argsort(pvals)
    sorted_p = pvals[sort_idx]
    qvals = np.zeros(n)
    for i in range(n): qvals[i] = sorted_p[i] * n / (i + 1)
    for i in range(n-2, -1, -1): qvals[i] = min(qvals[i], qvals[i+1])
    qvals = np.minimum(qvals, 1.0)
    res = np.zeros(n)
    res[sort_idx] = qvals
    return res

def significance_symbol(q):
    if q < 0.001: return "***"
    elif q < 0.01: return "**"
    elif q < 0.05: return "*"
    else: return "ns"

def annotate_stacked_percent(ax, pct_matrix, xpos, min_height=6.0):
    """Adds percentage text to stacked bars, hiding values below min_height to prevent overlapping."""
    for i in range(pct_matrix.shape[0]):
        bottom = 0
        for j, h in enumerate(pct_matrix[i]):
            if h >= min_height:
                ax.text(xpos[i], bottom + h/2, f"{h:.0f}%", ha="center", va="center", 
                        fontsize=11, fontweight="bold", color="white",
                        path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
            bottom += h

def annotate_bar_counts(ax, rects):
    """Automatically adds the exact count above each bar, rotated at 45 degrees to avoid overlap."""
    y_max = ax.get_ylim()[1]
    for rect in rects:
        h = rect.get_height()
        if h > 0:
            ax.text(rect.get_x() + rect.get_width()/2, h + (y_max * 0.03),
                    f"{int(h)}", ha="left", va="bottom", fontsize=10, fontweight="bold", rotation=45)

def _prettify_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def auto_generate_sites(layers_dir, out_csv):
    seqs, exposed_sets = {}, {}
    for iso in ISO_LABELS:
        pdb_path = os.path.join(layers_dir, f"{iso}.pdb")
        seq = []
        if os.path.exists(pdb_path):
            with open(pdb_path) as fh:
                for line in fh:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        rnum, rname = line[22:26].strip(), line[17:20].strip().upper()
                        if not seq or seq[-1][0] != rnum: seq.append((rnum, rname))
        seqs[iso] = seq
        _, _, exposed_sets[iso], _ = parse_pdb_stats(os.path.join(layers_dir, f"{iso}_surf1.pdb"))
        
    max_len = max(len(s) for s in seqs.values()) if seqs else 0
    alignment = {}
    for iso in ISO_LABELS:
        seq, offset = seqs[iso], max_len - len(seqs[iso])
        for idx, (rnum, rname) in enumerate(seq):
            aln_pos = idx + offset
            if aln_pos not in alignment: alignment[aln_pos] = {}
            is_exp, is_tit = (rnum in exposed_sets[iso]), (rname in TITRATABLE)
            alignment[aln_pos][iso] = {'Res': AA_3_TO_1.get(rname, "X"), 'Exposed': 'Yes' if is_exp else 'No', 'Titratable': 'Yes' if is_tit else 'No'}
            
    sites_rows = []
    for aln_pos in sorted(alignment.keys()):
        col_data = alignment[aln_pos]
        if not any(iso in col_data and col_data[iso]['Exposed'] == 'Yes' and col_data[iso]['Titratable'] == 'Yes' for iso in ISO_LABELS): continue
        res_set = set(col_data.get(iso, {}).get('Res', '-') for iso in ISO_LABELS)
        if len(res_set) > 1: 
            row = {"Aln_Pos": aln_pos}
            for iso in ISO_LABELS:
                if iso in col_data:
                    row.update({f"{iso}_Res": col_data[iso]['Res'], f"{iso}_Exposed": col_data[iso]['Exposed'], f"{iso}_Titratable": col_data[iso]['Titratable']})
                else:
                    row.update({f"{iso}_Res": "-", f"{iso}_Exposed": "No", f"{iso}_Titratable": "No"})
            sites_rows.append(row)
            
    if len(sites_rows) > 60:
        step = len(sites_rows) / 60.0
        sites_rows = [sites_rows[int(i*step)] for i in range(60)]
            
    if sites_rows:
        fieldnames = ["Aln_Pos"] + [f"{iso}_{t}" for iso in ISO_LABELS for t in ["Res", "Exposed", "Titratable"]]
        with open(out_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(sites_rows)
    return sites_rows


# ------------------------------------------------------------
# PLOTTING MODULE 1: COMPOSITION (4-PANEL)
# ------------------------------------------------------------
def plot_composition_4panel(labels, display_labels, atoms_prot, atoms_surf, residues_prot, residues_surf, 
                            frac_res, prot_pct, surf_pct, chi_gof_p, delta_charged, q_enrich, outdir):
    fig = plt.figure(figsize=(14, 16))
    gs  = GridSpec(3, 2, height_ratios=[1.2, 1.2, 0.30], hspace=1.0, wspace=0.45, figure=fig)

    axA = fig.add_subplot(gs[0,0])
    axB = fig.add_subplot(gs[0,1])
    axC = fig.add_subplot(gs[1,0])
    axD = fig.add_subplot(gs[1,1])
    axNote = fig.add_subplot(gs[2,:]); axNote.axis("off")
    for ax in [axA, axB, axC, axD]: _prettify_axes(ax)

    x = np.arange(len(labels))
    w = 0.18

    # Panel A
    rects1 = axA.bar(x - w*1.5, atoms_prot,    w, color=COLOR_PROTEIN,       label="Atoms (Protein)", edgecolor="black", linewidth=0.5)
    rects2 = axA.bar(x - w*0.5, atoms_surf,    w, color=COLOR_SURFACE,       label="Atoms (Surface)", edgecolor="black", linewidth=0.5)
    rects3 = axA.bar(x + w*0.5, residues_prot, w, color=COLOR_PROTEIN_LIGHT, label="Residues (Protein)", edgecolor="black", linewidth=0.5)
    rects4 = axA.bar(x + w*1.5, residues_surf, w, color=COLOR_SURFACE_LIGHT, label="Residues (Surface)", edgecolor="black", linewidth=0.5)
    axA.set_title("A) Total Atoms & Residues (Grouped)", pad=15, fontweight="bold", loc="left")
    axA.set_ylabel("Count")
    axA.set_xticks(x); axA.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    axA.set_ylim(0, max(atoms_prot) * 1.60 if atoms_prot else 100)
    annotate_bar_counts(axA, rects1); annotate_bar_counts(axA, rects2)
    annotate_bar_counts(axA, rects3); annotate_bar_counts(axA, rects4)
    axA.legend(frameon=False, ncols=2, loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=11)

    # Panel B
    bars = axB.bar(x, frac_res, width=0.5, color=COLOR_SURFACE, edgecolor="black", linewidth=0.5)
    axB.set_title("B) Fraction of Residues on Surface (%)", pad=15, fontweight="bold", loc="left")
    axB.set_ylabel("Percent"); axB.set_ylim(0, 120)
    axB.set_xticks(x); axB.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    for i, b in enumerate(bars):
        axB.text(b.get_x() + b.get_width()/2, b.get_height() + 2.0, f"{frac_res.iloc[i]:.1f}%", ha="center", va="bottom", fontsize=12, fontweight="bold")

    # Panel C
    bottom = np.zeros(len(labels))
    for j, c in enumerate(CLASS_LABELS):
        axC.bar(x, prot_pct[:, j], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="black", linewidth=0.5)
        bottom += prot_pct[:, j]
    axC.set_title("C) Composition — Protein (%)", pad=15, fontweight="bold", loc="left")
    axC.set_ylabel("Percent"); axC.set_ylim(0, 105)
    axC.set_xticks(x); axC.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    annotate_stacked_percent(axC, prot_pct, x)

    trans = btf(axC.transData, axC.transAxes)
    for i in range(len(labels)):
        axC.text(i, -0.30, f"χ² p={chi_gof_p[i]:.3f}", transform=trans, ha="center", va="top", fontsize=11, color="0.30",
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.80"), clip_on=False)

    # Panel D
    bottom = np.zeros(len(labels))
    for j, c in enumerate(CLASS_LABELS):
        axD.bar(x, surf_pct[:, j], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="black", linewidth=0.5)
        bottom += surf_pct[:, j]
    axD.set_title("D) Composition — Surface (%)", pad=15, fontweight="bold", loc="left")
    axD.set_ylabel("Percent"); axD.set_ylim(0, 120)
    axD.set_xticks(x); axD.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    annotate_stacked_percent(axD, surf_pct, x)

    for i in range(len(labels)):
        star  = significance_symbol(q_enrich[i])
        arrow = "↑" if delta_charged[i] > 0 else ("↓" if delta_charged[i] < 0 else "•")
        axD.text(i, 108, f"{arrow}{abs(delta_charged[i]):.1f}% {star}", ha="center", va="center", fontsize=12, color="black", fontweight="bold")

    bullets = [
        "Surface: LAYERS surf1.pdb; ATOM only; residues unique by (resname, chain, resid).",
        "Stacks show residue-class composition; numbers are segment percentages.",
        "Panel C: χ² goodness-of-fit (Surface vs Protein) p-values per isoform shown under bars.",
        "Panel D: Δcharged% = Surface – Protein; Fisher’s exact test with BH–FDR (stars)."
    ]
    axNote.text(0.5, 0.50, "• " + "\n• ".join(bullets), ha="center", va="center", fontsize=12,
                bbox=dict(boxstyle="round,pad=0.8", facecolor="#FAFAFA", edgecolor="#CFCFCF"))

    handles = [plt.Rectangle((0,0), 1, 1, color=CLASS_COLORS[c]) for c in CLASS_LABELS]
    fig.legend(handles, CLASS_LABELS, ncols=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.05), fontsize=13)

    fig.suptitle("Surface Composition of HK Isoforms (LAYERS Algorithm)", y=0.96, fontsize=18, fontweight="bold")

    out_png = os.path.join(outdir, "hk_surface_composition_4panel.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"-> Generated {out_png}")
    plt.close(fig)


# ------------------------------------------------------------
# PLOTTING MODULE 2: PROPENSITY & ASA (2-PANEL)
# ------------------------------------------------------------
def plot_propensity_asa_2panel(labels, display_labels, propensity_vals, asa_total, asa_tit, asa_perc, outdir):
    fig = plt.figure(figsize=(14, 7.5))
    gs = GridSpec(1, 2, width_ratios=[1.15, 1], wspace=0.35, figure=fig)

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    _prettify_axes(axA); _prettify_axes(axB)

    x = np.arange(len(labels))

    # Panel A: Propensity
    maxabs = max(abs(min(propensity_vals)), abs(max(propensity_vals))) if propensity_vals else 1.0
    if maxabs == 0: maxabs = 0.1 
    norm = TwoSlopeNorm(vmin=-maxabs-0.01, vcenter=0.0, vmax=maxabs+0.01)
    cmap = plt.get_cmap(COLOR_MAP_PROP)

    bars_prop = axA.bar(x, propensity_vals, color=cmap(norm(propensity_vals)), edgecolor="black", linewidth=0.8)
    axA.axhline(0, color="black", linewidth=1.0, linestyle="--")
    
    axA.set_ylim(-maxabs * 1.4, maxabs * 1.4)
    axA.set_xticks(x); axA.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    axA.set_ylabel("Log(Titratable on Surface / Total on Surface)")
    axA.set_title("A) Surface flip log-ratio per isoform\n(titratable side-chains)", loc="left", fontweight="bold", pad=15)
    axA.yaxis.grid(True, linestyle=":", linewidth=0.7, alpha=0.5)

    for i, rect in enumerate(bars_prop):
        v = propensity_vals[i]
        # FIXED: Enforce a firm 15% absolute padding so values near zero never overlap with the x-axis
        offset = 0.15 * maxabs if v >= 0 else -0.15 * maxabs
        axA.text(rect.get_x() + rect.get_width() / 2, v + offset, f"{v:.3f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=12, fontweight="bold")

    cax = inset_axes(axA, width="60%", height="4%", loc="lower center", bbox_to_anchor=(0.5, -0.28, 0.6, 0.12), bbox_transform=axA.transAxes)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Log propensity score", fontsize=12)

    # Panel B: ASA
    width = 0.36
    axB.bar(x - width / 2, asa_total, width=width, edgecolor="black", color=COLOR_PROTEIN, linewidth=0.8, label="Total Surface ASA (Å²)")
    axB.bar(x + width / 2, asa_tit, width=width, edgecolor="black", color=COLOR_SURFACE, linewidth=0.8, label="Titratable Surface ASA (Å²)")
    
    if asa_total and max(asa_total) > 0: axB.set_ylim(0, max(asa_total) * 1.25)
    axB.set_xticks(x); axB.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    
    for i, p in enumerate(asa_perc):
        if asa_total[i] > 0: axB.text(x[i] + width / 2, asa_tit[i] + (0.06 * max(asa_total)), f"{p:.1f}%", ha="center", fontsize=12, fontweight="bold")
    
    axB.set_ylabel("Accessible surface area (Å²)")
    axB.set_title("B) Accessible Surface Area (ASA)", loc="left", fontweight="bold", pad=15)
    axB.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), borderaxespad=0.0, ncol=1, fontsize=12)

    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.95]) 
    out_png = os.path.join(outdir, "hk_propensity_asa_2panel.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"-> Generated {out_png}")
    plt.close(fig)


# ------------------------------------------------------------
# PLOTTING MODULE 3: HEATMAP (1-PANEL)
# ------------------------------------------------------------
def plot_heatmap_1panel(labels, sites_rows, outdir):
    fig = plt.figure(figsize=(14.0, 9.5)) 
    gs  = GridSpec(2, 1, height_ratios=[1.0, 0.22], hspace=0.38, figure=fig)
    ax  = fig.add_subplot(gs[0, 0])
    ax_note = fig.add_subplot(gs[1, 0]); ax_note.axis("off")

    if not sites_rows:
        ax.text(0.5, 0.5, "Sites CSV data empty.", ha="center", va="center", fontsize=12, color="red")
    else:
        n_cols, n_rows = len(sites_rows), len(labels)
        ax.set_xlim(0, n_cols); ax.set_ylim(0, n_rows)
        ax.invert_yaxis(); ax.set_yticks(np.arange(n_rows) + 0.5)
        ax.set_yticklabels(labels, rotation=0, fontweight="bold") 

        step = max(1, int(np.ceil(n_cols / 22)))
        x_ticks = np.arange(0, n_cols, step) + 0.5
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([str(sites_rows[i].get("Aln_Pos", i)) for i in range(0, n_cols, step)], rotation=90)
        ax.set_title("Per-Isoform Residues at Variable Sites (filled = exposed + titratable)", fontweight="bold", pad=15)
        ax.set_xlabel("Aligned Residue Position (C-term anchored)")
        ax.set_ylabel("Isoform")

        for j, r in enumerate(sites_rows):
            for i, iso in enumerate(ISO_LABELS.keys()):
                aa  = get_csv_val(r, iso, "Res") or "-"
                exp = (str(get_csv_val(r, iso, "Exposed")).strip().lower() in ["yes", "true", "1", "y"])
                tit = (str(get_csv_val(r, iso, "Titratable")).strip().lower() in ["yes", "true", "1", "y"])
                
                grp = tit_group_letter(aa) if (exp and tit) else None
                face = HEATMAP_COLORS.get(grp, "white")

                rect = Rectangle((j, i), 1, 1, facecolor=face, edgecolor="0.80", linewidth=0.8)
                ax.add_patch(rect)

                text_color = "white" if face != "white" else "black"
                path_eff = [pe.withStroke(linewidth=2, foreground="black")] if face != "white" else None
                ax.text(j + 0.5, i + 0.5, aa if aa != "-" else "–", ha="center", va="center", color=text_color, fontsize=11, fontweight="bold", path_effects=path_eff)

        for j in range(n_cols + 1): ax.plot([j, j], [0, n_rows], color="0.92", lw=0.6, zorder=0)
        for i in range(n_rows + 1): ax.plot([0, n_cols], [i, i], color="0.92", lw=0.6, zorder=0)
        ax.set_frame_on(False)
        ax.set_aspect('auto') 

        swatches = [("basic (H/K/R)", HEATMAP_COLORS["basic"]), ("acidic (D/E)", HEATMAP_COLORS["acidic"]), ("other (Y/C)", HEATMAP_COLORS["other"])]
        for k, (lab, col) in enumerate(swatches):
            x0 = 0.05 + k * 0.25
            # FIXED: Pushed swatches fully below the X-axis for standalone heatmap
            ax.add_patch(Rectangle((x0, -0.45), 0.025, 0.12, transform=ax.transAxes, facecolor=col, edgecolor="black", lw=0.8, clip_on=False))
            ax.text(x0 + 0.035, -0.39, lab, transform=ax.transAxes, va="center", fontsize=12, fontweight="bold")

    ax_note.text(0.5, 0.5, "How to read: rows = isoforms; columns = variable alignment positions.\n"
                           "Filled color only when residue is both surface-exposed and titratable.\n"
                           "Colors: blue = basic (H/K/R); red = acidic (D/E); purple = Tyr/Cys.",
                 ha="center", va="center", fontsize=11, bbox=dict(boxstyle="round,pad=0.8", facecolor="#FAFAFA", edgecolor="#CFCFCF"))

    fig.suptitle("HK Isoforms — Surface-Exposed, Titratable, and Variable Sites", y=0.96, fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.95])
    
    out_png = os.path.join(outdir, "hk_variable_sites_heatmap.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"-> Generated {out_png}")
    plt.close(fig)


# ------------------------------------------------------------
# PLOTTING MODULE 4: THE COMBINED MEGA PLOT (7-PANEL)
# ------------------------------------------------------------
def plot_mega_7panel(labels, display_labels, atoms_prot, atoms_surf, residues_prot, residues_surf, 
                     frac_res, prot_pct, surf_pct, chi_gof_p, delta_charged, q_enrich, 
                     propensity_vals, asa_total, asa_tit, asa_perc, sites_rows, outdir):
                     
    fig = plt.figure(figsize=(24, 24)) 
    gs = GridSpec(4, 3, height_ratios=[1.2, 1.2, 0.9, 0.15], wspace=0.35, hspace=0.95, figure=fig)

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])
    axD = fig.add_subplot(gs[1, 0])
    axE = fig.add_subplot(gs[1, 1])
    axF = fig.add_subplot(gs[1, 2])
    axG = fig.add_subplot(gs[2, :])
    axNote = fig.add_subplot(gs[3, :]); axNote.axis("off")

    for ax in [axA, axB, axC, axD, axE, axF]: _prettify_axes(ax)

    x = np.arange(len(labels))
    w = 0.18

    # --- Panel A: Atoms & Residues
    rects1 = axA.bar(x - w*1.5, atoms_prot,    w, color=COLOR_PROTEIN,       label="Atoms (Protein)", edgecolor="black", linewidth=0.5)
    rects2 = axA.bar(x - w*0.5, atoms_surf,    w, color=COLOR_SURFACE,       label="Atoms (Surface)", edgecolor="black", linewidth=0.5)
    rects3 = axA.bar(x + w*0.5, residues_prot, w, color=COLOR_PROTEIN_LIGHT, label="Residues (Protein)", edgecolor="black", linewidth=0.5)
    rects4 = axA.bar(x + w*1.5, residues_surf, w, color=COLOR_SURFACE_LIGHT, label="Residues (Surface)", edgecolor="black", linewidth=0.5)
    axA.set_title("A) Total Atoms & Residues (Grouped)", loc="left", fontweight="bold", pad=10)
    axA.set_ylabel("Count")
    axA.set_xticks(x); axA.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    axA.set_ylim(0, max(atoms_prot) * 1.65 if atoms_prot else 100) 
    annotate_bar_counts(axA, rects1); annotate_bar_counts(axA, rects2)
    annotate_bar_counts(axA, rects3); annotate_bar_counts(axA, rects4)
    axA.legend(frameon=False, ncols=2, loc="upper center", bbox_to_anchor=(0.5, -0.28), fontsize=11)

    # --- Panel B: Fraction on Surface
    bars = axB.bar(x, frac_res, width=0.5, color=COLOR_SURFACE, edgecolor="black", linewidth=0.5)
    axB.set_title("B) Fraction of Residues on Surface (%)", loc="left", fontweight="bold", pad=10)
    axB.set_ylabel("Percent"); axB.set_ylim(0, 120)
    axB.set_xticks(x); axB.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    for i, b in enumerate(bars):
        axB.text(b.get_x() + b.get_width()/2, b.get_height() + 2.0, f"{frac_res.iloc[i]:.1f}%", ha="center", va="bottom", fontsize=12, fontweight="bold")

    # --- Panel C: Propensity
    maxabs = max(abs(min(propensity_vals)), abs(max(propensity_vals))) if propensity_vals else 1.0
    if maxabs == 0: maxabs = 0.1 
    norm = TwoSlopeNorm(vmin=-maxabs-0.01, vcenter=0.0, vmax=maxabs+0.01)
    cmap = plt.get_cmap(COLOR_MAP_PROP)
    bars_prop = axC.bar(x, propensity_vals, color=cmap(norm(propensity_vals)), edgecolor="black", linewidth=0.8)
    axC.axhline(0, color="black", linewidth=1.0, linestyle="--")
    axC.set_title("C) Surface flip log-ratio\n(titratable side-chains)", loc="left", fontweight="bold", pad=10)
    axC.set_ylabel("Log(Titratable Surf / Total Surf)")
    axC.set_xticks(x); axC.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    axC.set_ylim(-maxabs * 1.4, maxabs * 1.4) 

    for i, rect in enumerate(bars_prop):
        v = propensity_vals[i]
        # FIXED: Enforce a firm 15% absolute padding so values near zero never overlap with the x-axis line
        offset = 0.15 * maxabs if v >= 0 else -0.15 * maxabs
        axC.text(rect.get_x() + rect.get_width() / 2, v + offset, f"{v:.3f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=12, fontweight="bold")
    
    cax = inset_axes(axC, width="60%", height="4%", loc="lower center", bbox_to_anchor=(0.5, -0.32, 0.6, 0.08), bbox_transform=axC.transAxes)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Log propensity score", fontsize=11)

    # --- Panel D: Composition Protein
    bottom = np.zeros(len(labels))
    for j, c in enumerate(CLASS_LABELS):
        axD.bar(x, prot_pct[:, j], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="black", linewidth=0.5)
        bottom += prot_pct[:, j]
    axD.set_title("D) Composition — Protein (%)", loc="left", fontweight="bold", pad=10)
    axD.set_ylabel("Percent"); axD.set_ylim(0, 105)
    axD.set_xticks(x); axD.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    annotate_stacked_percent(axD, prot_pct, x)
    trans = btf(axD.transData, axD.transAxes)
    for i in range(len(labels)):
        axD.text(i, -0.38, f"χ² p={chi_gof_p[i]:.3f}", transform=trans, ha="center", va="top", fontsize=11, color="0.30",
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.80"), clip_on=False)

    # --- Panel E: Composition Surface
    bottom = np.zeros(len(labels))
    for j, c in enumerate(CLASS_LABELS):
        axE.bar(x, surf_pct[:, j], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="black", linewidth=0.5)
        bottom += surf_pct[:, j]
    axE.set_title("E) Composition — Surface (%)", loc="left", fontweight="bold", pad=10)
    axE.set_ylabel("Percent"); axE.set_ylim(0, 120)
    axE.set_xticks(x); axE.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    annotate_stacked_percent(axE, surf_pct, x)
    for i in range(len(labels)):
        star  = significance_symbol(q_enrich[i])
        arrow = "↑" if delta_charged[i] > 0 else ("↓" if delta_charged[i] < 0 else "•")
        axE.text(i, 108, f"{arrow}{abs(delta_charged[i]):.1f}% {star}", ha="center", va="center", fontsize=12, color="black", fontweight="bold")

    # --- Panel F: ASA Bars
    width = 0.36
    axF.bar(x - width / 2, asa_total, width=width, edgecolor="black", color=COLOR_PROTEIN, linewidth=0.8, label="Total Surface ASA (Å²)")
    axF.bar(x + width / 2, asa_tit, width=width, edgecolor="black", color=COLOR_SURFACE, linewidth=0.8, label="Titratable Surface ASA (Å²)")
    axF.set_title("F) Accessible Surface Area (ASA)", loc="left", fontweight="bold", pad=10)
    axF.set_ylabel("Accessible surface area (Å²)")
    if max(asa_total) > 0: axF.set_ylim(0, max(asa_total) * 1.25)
    axF.set_xticks(x); axF.set_xticklabels(display_labels, rotation=0, ha="center", fontweight="bold")
    for i, p in enumerate(asa_perc):
        if asa_total[i] > 0: axF.text(x[i] + width / 2, asa_tit[i] + (0.06 * max(asa_total)), f"{p:.1f}%", ha="center", fontsize=12, fontweight="bold")
    axF.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), borderaxespad=0.0, ncol=1, fontsize=11)

    # --- Panel G: Variable Sites Heatmap
    if sites_rows:
        n_cols, n_rows = len(sites_rows), len(ISO_LABELS)
        axG.set_xlim(0, n_cols); axG.set_ylim(0, n_rows)
        axG.invert_yaxis(); axG.set_yticks(np.arange(n_rows) + 0.5)
        axG.set_yticklabels(labels, rotation=0, fontweight="bold") 
        
        axG.set_aspect('auto') 

        step = max(1, int(np.ceil(n_cols / 18)))
        x_ticks = np.arange(0, n_cols, step) + 0.5
        axG.set_xticks(x_ticks)
        axG.set_xticklabels([str(sites_rows[i].get("Aln_Pos", i)) for i in range(0, n_cols, step)], rotation=90)
        axG.set_title("G) Per-Isoform Residues at Variable Sites (filled = exposed + titratable)", loc="left", fontweight="bold", pad=10)
        axG.set_xlabel("Aligned Residue Position (C-term anchored)")

        for j, r in enumerate(sites_rows):
            for i, iso in enumerate(ISO_LABELS.keys()):
                aa  = get_csv_val(r, iso, "Res") or "-"
                exp = (str(get_csv_val(r, iso, "Exposed")).strip().lower() in ["yes", "true", "1", "y"])
                tit = (str(get_csv_val(r, iso, "Titratable")).strip().lower() in ["yes", "true", "1", "y"])
                grp = tit_group_letter(aa) if (exp and tit) else None
                face = HEATMAP_COLORS.get(grp, "white")
                axG.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="0.80", linewidth=0.8))
                
                text_color = "white" if face != "white" else "black"
                path_eff = [pe.withStroke(linewidth=1.5, foreground="black")] if face != "white" else None
                axG.text(j + 0.5, i + 0.5, aa if aa != "-" else "–", ha="center", va="center", color=text_color, fontsize=9, fontweight="bold", path_effects=path_eff)

        for j in range(n_cols + 1): axG.plot([j, j], [0, n_rows], color="0.92", lw=0.6, zorder=0)
        for i in range(n_rows + 1): axG.plot([0, n_cols], [i, i], color="0.92", lw=0.6, zorder=0)
        axG.set_frame_on(False)
        
        swatches = [("basic (H/K/R)", HEATMAP_COLORS["basic"]), ("acidic (D/E)", HEATMAP_COLORS["acidic"]), ("other (Y/C)", HEATMAP_COLORS["other"])]
        for k, (lab, col) in enumerate(swatches):
            x0 = 0.05 + k * 0.25
            # FIXED: Pushed heavily down (-0.65) to completely clear the 90-degree rotated X-axis labels
            axG.add_patch(Rectangle((x0, -0.65), 0.025, 0.15, transform=axG.transAxes, facecolor=col, edgecolor="black", lw=0.8, clip_on=False))
            axG.text(x0 + 0.035, -0.57, lab, transform=axG.transAxes, va="center", fontsize=12, fontweight="bold")

    # --- Note Box & Master Legends
    bullets = [
        "Panel A & B: Total atoms/residues and fraction exposed on the surface (LAYERS algorithm).",
        "Panel C: Surface flip log-ratio (propensity) of titratable side-chains (values > 0 indicate enrichment).",
        "Panel D & E: Residue-class composition. χ² goodness-of-fit p-values shown under D. Δcharged% with Fisher’s exact test (BH–FDR stars) shown above E.",
        "Panel F: Total vs. Titratable Accessible Surface Area (ASA) calculated via NACCESS.",
        "Panel G: Categorical sequence alignment of variable sites. Highlighted cells indicate exposed and titratable residues."
    ]
    axNote.text(0.5, 0.70, "• " + "\n• ".join(bullets), ha="center", va="center", fontsize=13, bbox=dict(boxstyle="round,pad=1.0", facecolor="#FAFAFA", edgecolor="#CFCFCF"))

    handles = [plt.Rectangle((0,0), 1, 1, color=CLASS_COLORS[c]) for c in CLASS_LABELS]
    fig.legend(handles, CLASS_LABELS, ncols=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.02), fontsize=13)

    fig.suptitle("Comprehensive Surface Composition & Variable Site Analysis of HK Isoforms", y=0.96, fontsize=20, fontweight="bold")
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
    
    out_png = os.path.join(outdir, "hk_comprehensive_mega_plot.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"-> Generated {out_png}")
    plt.close(fig)


# ------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--layers", default=".", help="Directory containing *_surf1.pdb and .pdb files")
    p.add_argument("--rsa", default=".", help="Directory containing .rsa files")
    p.add_argument("--sites", default="hk_layers_sites.csv", help="Auto-generated CSV filename")
    p.add_argument("--outdir", default="hk_plots", help="Output directory")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    
    labels, display_labels = [], []
    atoms_prot, atoms_surf, residues_prot, residues_surf = [], [], [], []
    prot_pct_list, surf_pct_list, chi_gof_p, delta_charged, fisher_p_raw = [], [], [], [], []
    propensity_vals, asa_total, asa_tit, asa_perc = [], [], [], []

    print("Analyzing Hexokinase structures...")
    for iso_key, iso_label in ISO_LABELS.items():
        labels.append(iso_label)
        display_labels.append(iso_label.replace(" (", "\n(")) # Two neat lines
        
        full_pdb = os.path.join(args.layers, f"{iso_key}.pdb")
        surf_pdb = os.path.join(args.layers, f"{iso_key}_surf1.pdb")
        
        a_prot, r_prot_set, p_resnum_set, p_resnum_map = parse_pdb_stats(full_pdb)
        a_surf, r_surf_set, s_resnum_set, _ = parse_pdb_stats(surf_pdb)
        
        atoms_prot.append(a_prot); atoms_surf.append(a_surf)
        residues_prot.append(len(r_prot_set)); residues_surf.append(len(r_surf_set))
        
        p_pcts, p_counts = get_composition(r_prot_set)
        s_pcts, s_counts = get_composition(r_surf_set)
        prot_pct_list.append(p_pcts); surf_pct_list.append(s_pcts)
        delta_charged.append((s_pcts[2] + s_pcts[3]) - (p_pcts[2] + p_pcts[3]))
        
        if sum(s_counts) > 0 and sum(p_counts) > 0:
            f_obs = s_counts
            f_exp = np.where(np.array(p_counts) == 0, 1e-6, np.array(p_counts) / sum(p_counts) * sum(s_counts))
            _, chi_p = chisquare(f_obs, f_exp)
            chi_gof_p.append(chi_p)
            
            surf_chg = s_counts[2] + s_counts[3]; surf_unchg = s_counts[0] + s_counts[1]
            int_chg = max(0, p_counts[2] + p_counts[3] - surf_chg); int_unchg = max(0, p_counts[0] + p_counts[1] - surf_unchg)
            _, fish_p = fisher_exact([[surf_chg, surf_unchg], [int_chg, int_unchg]])
            fisher_p_raw.append(fish_p)
        else:
            chi_gof_p.append(1.0); fisher_p_raw.append(1.0)
            
        flips_on_surf = sum(1 for k in s_resnum_set if (k in p_resnum_map and p_resnum_map[k] in TITRATABLE))
        frac_surf = (flips_on_surf + PSEUDOCOUNT) / (len(s_resnum_set) + PSEUDOCOUNT) if len(s_resnum_set) > 0 else PSEUDOCOUNT
        propensity_vals.append(np.log(frac_surf))
        
        rsa_file = os.path.join(args.rsa, f"{iso_key}.rsa")
        if not os.path.exists(rsa_file): rsa_file = os.path.join(args.rsa, f"{iso_key.upper()}.rsa")
        rsa_map = parse_rsa(rsa_file)
        
        t_asa, t_tit = 0.0, 0.0
        for resnum in s_resnum_set:
            if resnum in rsa_map:
                resname, abs_val = rsa_map[resnum]
                t_asa += abs_val
                if resname in TITRATABLE: t_tit += abs_val
        asa_total.append(t_asa); asa_tit.append(t_tit)
        asa_perc.append((t_tit / t_asa * 100.0) if t_asa > 0 else 0.0)

    q_enrich = bh_fdr(fisher_p_raw)
    frac_res = pd.Series([(rs/rp*100) if rp > 0 else 0 for rs, rp in zip(residues_surf, residues_prot)], index=labels)
    prot_pct, surf_pct = np.array(prot_pct_list), np.array(surf_pct_list)

    sites_rows = auto_generate_sites(args.layers, args.sites)

    print("\nGenerating pristine, color-blind friendly plots...")
    
    # Generate the 3 separated plots
    plot_composition_4panel(labels, display_labels, atoms_prot, atoms_surf, residues_prot, residues_surf, frac_res, prot_pct, surf_pct, chi_gof_p, delta_charged, q_enrich, args.outdir)
    plot_propensity_asa_2panel(labels, display_labels, propensity_vals, asa_total, asa_tit, asa_perc, args.outdir)
    plot_heatmap_1panel(labels, sites_rows, args.outdir)
    
    # Generate the combined mega-plot
    plot_mega_7panel(labels, display_labels, atoms_prot, atoms_surf, residues_prot, residues_surf, frac_res, prot_pct, surf_pct, chi_gof_p, delta_charged, q_enrich, propensity_vals, asa_total, asa_tit, asa_perc, sites_rows, args.outdir)
    
    print("\nALL DONE! Check your output directory for publication-ready figures.")

if __name__ == "__main__":
    main()
