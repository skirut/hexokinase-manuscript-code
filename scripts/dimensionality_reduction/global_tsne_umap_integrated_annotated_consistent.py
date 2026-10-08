"""Global t-SNE and UMAP analysis of Hexokinase RMSD trajectories.

Purpose
-------
Builds a combined dataset from apo and G6P-bound RMSD trajectories for
HK1-HK4 and the modeled mutants at pH 3.0, 5.0, and 7.0. Each RMSD
trajectory is interpolated to 1000 points, followed by t-SNE and UMAP
projection and a publication-style comparison plot.

Inputs
------
- ``$HEXOKINASE_DATA_ROOT/structures/<system>_pH<pH>/prod/rmsd_analysis.dat``
- ``$HEXOKINASE_DATA_ROOT/G6P_structure/<bound-system>_pH<pH>/prod/rmsd_analysis.dat``
- Bound-system mappings are defined in ``system_mapping``.

Outputs
-------
- ``results/dimensionality_reduction/Final_Comparison_TSNE_UMAP_CB.png``

Run from
--------
Run from the repository root with ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/dimensionality_reduction/global_tsne_umap_integrated_annotated_consistent.py``

Dependencies
------------
Python, NumPy, pandas, Matplotlib, scikit-learn, umap-learn, optional
adjustText, and ``scripts.utils.paths``.
"""

import os
import io
import string
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from sklearn.manifold import TSNE
import umap.umap_ as umap

# Try to import adjustText for better annotation placement
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

warnings.filterwarnings("ignore")

# ==========================================
# 1. DIRECTORY SETUP & MAPPINGS
# ==========================================
from scripts.utils.paths import APO_DIR, BOUND_DIR, DIM_REDUCTION_DIR

apo_dir = APO_DIR
bound_dir = BOUND_DIR
output_dir = DIM_REDUCTION_DIR
os.makedirs(output_dir, exist_ok=True)

display_names = {
    '1qha': 'HK1', '1qha_G87E': 'HK1_G87E',
    '2nzt': 'HK2', '2nzt_D84Y': 'HK2_D84Y', '2nzt_D861Y': 'HK2_D861Y', '2nzt_T680M': 'HK2_T680M',
    '3hm8': 'HK3', '3hm8_D867N': 'HK3_D867N', '3hm8_D867Y': 'HK3_D867Y', '3hm8_T546M': 'HK3_T546M',
    '1v4s': 'HK4', '1v4s_D205N': 'HK4_D205N'
}

system_mapping = {
    '1qha': '1QHA_model_with_G6P', '1qha_G87E': '1QHA_G87E_G6P',
    '2nzt': '2nzt_model_with_G6P', '2nzt_D84Y': '2nzt_D84Y_G6P',
    '2nzt_D861Y': '2nzt_D861Y_G6P', '2nzt_T680M': '2nzt_T680M_G6P',
    '3hm8': '3hm8_model_with_G6P', '3hm8_D867N': '3hm8_D867N_G6P',
    '3hm8_D867Y': '3hm8_D867Y_G6P', '3hm8_T546M': '3hm8_T546M_G6P',
    '1v4s': None, '1v4s_D205N': None
}

pH_levels = ['3.0', '5.0', '7.0']

# ==========================================
# 2. COLOR-BLIND FRIENDLY AESTHETICS (Wong Palette)
# ==========================================
# Optimized for Deuteranopia, Protanopia, and Tritanopia
CB_COLORS = {
    'HK1': '#E69F00', # Orange
    'HK2': '#56B4E9', # Sky Blue
    'HK3': '#009E73', # Bluish Green
    'HK4': '#CC79A7'  # Reddish Purple
}

MARKER_MAP = {'HK1': 'o', 'HK2': 's', 'HK3': '^', 'HK4': 'D'}

# ==========================================
# 3. DATA EXTRACTION (Simplified for RMSD)
# ==========================================
def extract_data(filepath):
    if not os.path.exists(filepath): return None
    try:
        df = pd.read_csv(filepath, sep=r'\s+', comment='#', header=None)
        y = df[1].values if df.shape[1] > 1 else df[0].values
        return np.interp(np.linspace(0, 1, 1000), np.linspace(0, 1, len(y)), y)
    except: return None

def compile_dataset():
    data, labels = [], []
    for sys_int, sys_disp in display_names.items():
        bound_prefix = system_mapping[sys_int]
        for ph in pH_levels:
            for state in ['Apo', 'Bound']:
                if state == 'Bound' and not bound_prefix: continue
                path = os.path.join(bound_dir if state == 'Bound' else apo_dir, 
                                    f"{bound_prefix if state == 'Bound' else sys_int}_pH{ph}", 
                                    'prod', 'rmsd_analysis.dat')
                y = extract_data(path)
                if y is not None:
                    data.append(y)
                    labels.append({
                        'Enzyme': sys_disp.split('_')[0],
                        'Full_Name': sys_disp,
                        'pH': ph,
                        'State': state,
                        'ID': f"{sys_disp} {state} pH{ph}"
                    })
    return np.array(data), pd.DataFrame(labels)

# ==========================================
# 4. PLOTTING SINGLE FIGURE (t-SNE vs UMAP)
# ==========================================
def plot_comparison(embeddings, df_labels):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    methods = [('t-SNE', 'a'), ('UMAP', 'b')]

    for ax, (method, letter) in zip(axes, methods):
        X_2d = embeddings[method]
        texts = []
        
        for i in range(len(df_labels)):
            enz = df_labels['Enzyme'].iloc[i]
            state = df_labels['State'].iloc[i]
            color = CB_COLORS[enz]
            marker = MARKER_MAP[enz]
            
            # Scatter plot
            ax.scatter(X_2d[i, 0], X_2d[i, 1], s=180, marker=marker,
                       facecolors=color if state == 'Bound' else 'white',
                       edgecolors=color, linewidths=2.5, alpha=0.8)
            
            # TARGETED ANNOTATION: Only HK1 Bound (WT)
            if df_labels['Full_Name'].iloc[i] == 'HK1' and state == 'Bound':
                label_str = f"HK1 Bnd pH{df_labels['pH'].iloc[i]}"
                texts.append(ax.text(X_2d[i, 0], X_2d[i, 1], label_str, 
                                     fontsize=10, fontweight='bold', color='black'))

        if HAS_ADJUST_TEXT:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=1))

        # Aesthetics
        ax.set_title(f"{method} Projection (RMSD)", fontsize=16, fontweight='bold')
        ax.set_xlabel(f"{method} 1", fontsize=12)
        ax.set_ylabel(f"{method} 2", fontsize=12)
        ax.text(-0.05, 1.05, letter, transform=ax.transAxes, fontsize=24, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.3)

    # Global Legend
    legend_elements = [
        mlines.Line2D([], [], color='none', label='STATE:'),
        mlines.Line2D([], [], color='white', marker='o', mfc='gray', mec='black', label='Bound'),
        mlines.Line2D([], [], color='white', marker='o', mfc='white', mec='gray', mew=2, label='Apo'),
        mlines.Line2D([], [], color='none', label=' '),
        mlines.Line2D([], [], color='none', label='ISOFORMS:')
    ]
    for enz, col in CB_COLORS.items():
        legend_elements.append(mlines.Line2D([], [], color='white', marker=MARKER_MAP[enz], 
                                             mfc=col, mec='black', label=enz))

    axes[1].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
    
    plt.tight_layout()
    out_path = os.path.join(output_dir, "Final_Comparison_TSNE_UMAP_CB.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Success! Figure saved to: {out_path}")

# ==========================================
# 5. EXECUTION
# ==========================================
if __name__ == "__main__":
    X, df = compile_dataset()
    print(f"Data Compiled: {len(df)} trajectories.")
    
    # Compute Projections
    tsne_res = TSNE(n_components=2, perplexity=min(30, len(df)//2), random_state=42).fit_transform(X)
    umap_res = umap.UMAP(n_components=2, n_neighbors=min(15, len(df)-1), random_state=42).fit_transform(X)
    
    plot_comparison({'t-SNE': tsne_res, 'UMAP': umap_res}, df)

