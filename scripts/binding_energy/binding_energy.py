"""Hexokinase binding-energy master-grid plotting.

Purpose
-------
Reads VMD-generated ``binding_energy.dat`` files for the bound
Hexokinase systems and summarizes the total interaction energy as
Elec + VdW. The script organizes HK1-HK4 wild-type and mutant systems
across pH 3.0, 5.0, and 7.0 and produces one publication-ready grid.

Inputs
------
- ``$HEXOKINASE_DATA_ROOT/G6P_structure/*_pH*/prod/binding_energy.dat``
- Folder names must match the HK_MAPPING defined below.
- The input table must contain Elec, VdW, and Time columns (case
  normalization is handled by the script).

Outputs
-------
- ``results/binding_energy/Figure_Binding_Energy_Master_Grid.png``

Run from
--------
Run from the repository root so that ``scripts.utils.paths`` can be
imported. ``HEXOKINASE_DATA_ROOT`` must point to the study data root.

Usage
-----
``python scripts/binding_energy/binding_energy.py``

Dependencies
------------
Python, NumPy, pandas, Matplotlib, PyYAML, and this repository's
``scripts.utils.paths`` module.
"""

import os
import glob
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# CONFIG
# =========================
from scripts.utils.paths import BOUND_DIR, BINDING_ENERGY_DIR

BASE_DIR = BOUND_DIR
OUTDIR = BINDING_ENERGY_DIR
os.makedirs(OUTDIR, exist_ok=True)

# Time units and limits
TIME_UNITS = "ns"   
X_MAX_NS   = 300.0  

# Smoothed vs raw styling
RAW_ALPHA     = 0.35
RAW_LINEWIDTH = 1.0
MA_LINEWIDTH  = 2.0

# Global font scaling
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 12,
})

# Dictionary to reliably parse the Hexokinase folder structures
HK_MAPPING = {
    '1QHA_model_with_G6P': ('HK1', 'WT'),
    '1QHA_G87E_G6P':       ('HK1', 'G87E'),
    '2nzt_model_with_G6P': ('HK2', 'WT'),
    '2nzt_D84Y_G6P':       ('HK2', 'D84Y'),
    '2nzt_D861Y_G6P':      ('HK2', 'D861Y'),
    '2nzt_T680M_G6P':      ('HK2', 'T680M'),
    '3hm8_model_with_G6P': ('HK3', 'WT'),
    '3hm8_D867N_G6P':      ('HK3', 'D867N'),
    '3hm8_D867Y_G6P':      ('HK3', 'D867Y'),
    '3hm8_T546M_G6P':      ('HK3', 'T546M'),
    '1v4s_model_with_G6P': ('HK4', 'WT'),
    '1v4s_D205N_G6P':      ('HK4', 'D205N')
}

# Define strict grid order
VARIANT_ORDER = [
    'HK1_WT', 'HK1_G87E',
    'HK2_WT', 'HK2_D84Y', 'HK2_D861Y', 'HK2_T680M',
    'HK3_WT', 'HK3_D867N', 'HK3_D867Y', 'HK3_T546M',
    'HK4_WT', 'HK4_D205N'
]
PH_LEVELS = ['3.0', '5.0', '7.0']

# =========================
# HELPERS
# =========================
def parse_folder_meta(parent_name):
    if "_pH" not in parent_name:
        raise ValueError(f"Cannot parse pH from folder: {parent_name}")
        
    prefix, ph_token = parent_name.rsplit("_pH", 1)
    ph_value = ph_token.strip()
    
    base, mut = HK_MAPPING.get(prefix, (prefix, "WT"))
    variant = f"{base}_{mut}"
    
    return base, mut, ph_value, variant

def time_to_ns(series: pd.Series) -> pd.Series:
    return series / 1000.0 if TIME_UNITS.lower() == "ps" else series

def get_subplot_label(index):
    if index < 26:
        return string.ascii_uppercase[index]
    else:
        return string.ascii_uppercase[(index // 26) - 1] + string.ascii_uppercase[index % 26]

def collect_data():
    energy_files = glob.glob(os.path.join(BASE_DIR, "*_pH*", "prod", "binding_energy.dat"))
    print(f"Found {len(energy_files)} binding_energy.dat files")

    master_data = {v: {ph: None for ph in PH_LEVELS} for v in VARIANT_ORDER}

    for filepath in energy_files:
        try:
            parent = os.path.basename(os.path.dirname(os.path.dirname(filepath)))
            base, mut, ph_value, variant = parse_folder_meta(parent)
        except Exception as e:
            continue

        if variant not in master_data or ph_value not in PH_LEVELS:
            continue

        df = pd.read_csv(filepath, sep=r"\s+", comment="#", engine="python")
        
        cols_upper = [c.upper() for c in df.columns]
        if 'ELECT' in cols_upper and 'VDW' in cols_upper and 'TS' in cols_upper:
            df = df.rename(columns=dict(zip(df.columns, cols_upper)))
            df = df.rename(columns={'ELECT': 'Elec', 'VDW': 'VdW', 'TS': 'Time'})
            
        if not {"Elec", "VdW", "Time"}.issubset(df.columns):
            continue

        if "Frame" in df.columns:
            df = df[df["Frame"] != 0]

        df["Total"]    = df["Elec"] + df["VdW"]
        df["Time_ns"]  = time_to_ns(df["Time"])
        df["Total_MA"] = df["Total"].rolling(window=10, min_periods=1).mean()

        master_data[variant][ph_value] = df

    return master_data

# =========================
# MAIN PLOTTING
# =========================
def main():
    data = collect_data()
    
    # Filter out completely empty rows (variants with zero data across all pH levels)
    valid_variants = []
    for var in VARIANT_ORDER:
        if any(data[var][ph] is not None for ph in PH_LEVELS):
            valid_variants.append(var)
            
    n_rows = len(valid_variants)
    n_cols = len(PH_LEVELS)
    
    if n_rows == 0:
        print("No valid data found to plot.")
        return

    # Create the unified figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 3 * n_rows), sharex=True)
    plt.subplots_adjust(hspace=0.4, wspace=0.25, top=0.94, bottom=0.08)

    fig.suptitle("Time Evolution of Total Binding Energies", fontsize=24, fontweight='bold', y=0.98)
    
    plot_counter = 0

    for row_idx, variant in enumerate(valid_variants):
        # Calculate Y-limits for the current row to keep scale uniform across pHs for the same variant
        row_vals = []
        for ph in PH_LEVELS:
            df = data[variant][ph]
            if df is not None and not df.empty:
                row_vals.extend(df["Total"].tolist())
                
        if row_vals:
            y_min, y_max = min(row_vals), max(row_vals)
            pad = 0.1 * (y_max - y_min if y_max > y_min else 1.0)
            row_ylim = (y_min - pad, y_max + pad)
        else:
            row_ylim = (-100, 100) # Fallback

        for col_idx, ph in enumerate(PH_LEVELS):
            ax = axes[row_idx, col_idx] if n_rows > 1 else axes[col_idx]
            df = data[variant][ph]
            
            # Add Panel Label (A, B, C...)
            label_text = get_subplot_label(plot_counter)
            ax.text(-0.1, 1.05, label_text, transform=ax.transAxes, 
                    fontsize=16, fontweight='bold', va='bottom', ha='left')
            plot_counter += 1

            ax.set_title(f"{variant.replace('_', ' ')} | pH {ph}", fontsize=14, fontweight='bold')

            if df is not None and not df.empty:
                ax.plot(df["Time_ns"], df["Total"], alpha=RAW_ALPHA, color="#1f77b4", linewidth=RAW_LINEWIDTH, label="Raw Energy")
                ax.plot(df["Time_ns"], df["Total_MA"], color="#ff7f0e", linewidth=MA_LINEWIDTH, label="Running Average")
                
                ax.set_xlim(0, X_MAX_NS)
                ax.set_ylim(row_ylim)
                ax.grid(True, linestyle='--', alpha=0.5)

                # Mean and SD box
                avg = df["Total"].mean()
                std = df["Total"].std(ddof=1)
                ax.text(0.03, 0.05, f"{avg:.1f} ± {std:.1f} kcal/mol",
                        transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=12, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"))
            else:
                ax.text(0.5, 0.5, "Data Not Available", ha='center', va='center', 
                        transform=ax.transAxes, fontsize=12, color='gray')
                ax.set_xlim(0, X_MAX_NS)
                ax.set_ylim(-100, 100)
                ax.set_xticks([])
                ax.set_yticks([])

            # Axis labels
            if row_idx == n_rows - 1:
                ax.set_xlabel("Time (ns)", fontsize=13)
            if col_idx == 0:
                ax.set_ylabel("Binding Energy\n(kcal/mol)", fontsize=13)

            # Add legend to only the very first plot
            if row_idx == 0 and col_idx == 0 and df is not None:
                handles, labels = ax.get_legend_handles_labels()
                fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.95), frameon=False, fontsize=14)

    # Unified Explanatory Box at the bottom
    desc_text = (
        "Mean ± SD: Average and variability of Total binding energy (kcal/mol)\n"
        "Lower (more negative) values → stronger binding | Higher (less negative/positive) values → weaker binding"
    )
    fig.text(0.5, 0.02, desc_text, ha="center", va="center", fontsize=13, fontweight="medium",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", alpha=1.0, edgecolor="#cccccc"))

    outpath = os.path.join(OUTDIR, "Figure_Binding_Energy_Master_Grid.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[OK] Publication-ready master grid saved to {outpath}")

if __name__ == "__main__":
    main()

