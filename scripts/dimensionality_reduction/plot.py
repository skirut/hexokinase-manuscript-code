#!/usr/bin/env python3
"""Publication plot for the global Hexokinase conformational landscape.

Purpose
-------
Combines precomputed t-SNE and UMAP coordinates into a 2x2 figure
covering structural features and integrated structural-ligand features,
with state, isoform, and pH encoded in the plot.

Inputs
------
- ``Structural_4_Metrics_Coordinates.csv``
- ``Omni_8_Metrics_Coordinates.csv``
- The CSV files must be in the current working directory and contain
  the coordinate and annotation columns referenced by the script.

Outputs
-------
- ``Figure2_Global_Conformational_Landscape.png``
- ``Figure2_Global_Conformational_Landscape.tiff``
  Both are written to the current working directory.

Run from
--------
Run from the directory containing the two CSV input files. This script
does not use ``scripts.utils.paths`` for its input/output locations.

Usage
-----
``python scripts/dimensionality_reduction/plot.py``

Dependencies
------------
Python, pandas, Matplotlib, and optional adjustText.
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import matplotlib

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

# =========================
# STYLES (Okabe-Ito Palette)
# =========================
PH_MARKERS = {3.0: "o", 5.0: "^", 7.0: "s"}

ENZ_COLORS = {
    "HK1": "#D55E00",
    "HK2": "#0072B2",
    "HK3": "#009E73",
    "HK4": "#F0E442"
}

PH_COLORS = {
    3.0: "#E69F00",
    5.0: "#56B4E9",
    7.0: "#CC79A7"
}

# =========================
# MASTER LEGEND (Two-Row Layout)
# =========================
def add_master_legend(fig):

    MARKER_SIZE  = 18
    LEGEND_FONT  = 19
    HEADER_FONT  = 20

    # ---------- ROW 1 : STATE + ISOFORM ----------
    row1 = [
        Line2D([0], [0], color='none', label='STATE:'),

        Line2D(
            [0], [0],
            marker='s', color='w',
            markerfacecolor='gray',
            markeredgecolor='k',
            markersize=MARKER_SIZE,
            label='Bound (G6P)'
        ),

        Line2D(
            [0], [0],
            marker='s', color='none',
            markerfacecolor='none',
            markeredgecolor='gray',
            markeredgewidth=3.0,
            markersize=MARKER_SIZE,
            label='Unbound (Apo)'
        ),

        Line2D([0], [0], color='none', label='   |   '),

        Line2D([0], [0], color='none', label='ISOFORM:')
    ]

    for enz, col in ENZ_COLORS.items():
        row1.append(
            Line2D(
                [0], [0],
                marker='o', color='w',
                markerfacecolor=col,
                markeredgecolor='k',
                markersize=MARKER_SIZE,
                label=enz
            )
        )

    leg1 = fig.legend(
        handles=row1,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.065),
        ncol=len(row1),
        frameon=False,
        fontsize=LEGEND_FONT,
        handletextpad=0.5,
        columnspacing=1.2,
        borderaxespad=0.0
    )

    # ---------- ROW 2 : pH CONDITION ----------
    row2 = [Line2D([0], [0], color='none', label='pH CONDITION:')]

    for ph, mark in PH_MARKERS.items():
        row2.append(
            Line2D(
                [0], [0],
                marker=mark, color='w',
                markerfacecolor=PH_COLORS[ph],
                markeredgecolor='k',
                markersize=MARKER_SIZE,
                label=f"pH {ph}"
            )
        )

    leg2 = fig.legend(
        handles=row2,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.015),
        ncol=len(row2),
        frameon=False,
        fontsize=LEGEND_FONT,
        handletextpad=0.5,
        columnspacing=1.2,
        borderaxespad=0.0
    )

    # Bold section headers
    header_labels = {'STATE:', 'ISOFORM:', 'pH CONDITION:'}
    for leg in (leg1, leg2):
        for txt in leg.get_texts():
            if txt.get_text().strip() in header_labels:
                txt.set_fontsize(HEADER_FONT)
                txt.set_fontweight('bold')

    fig.add_artist(leg1)
    fig.add_artist(leg2)


# =========================
# PANEL PLOTTING
# =========================
def plot_panel(ax, df, x_col, y_col, title, panel_label):

    texts = []

    for _, row in df.iterrows():

        is_bound = row["State"] == "Bound"

        face_color = (
            ENZ_COLORS.get(row["Isoform"], "gray")
            if is_bound else "none"
        )

        edge_color = (
            "black"
            if is_bound else ENZ_COLORS.get(row["Isoform"], "gray")
        )

        lw = 1.0 if is_bound else 2.5

        ax.scatter(
            row[x_col],
            row[y_col],
            marker=PH_MARKERS.get(row["pH"], "o"),
            s=280,
            facecolors=face_color,
            edgecolors=edge_color,
            linewidths=lw,
            zorder=3
        )

        lbl = f"{row['Isoform']}_{row['State'][:3]}_pH{row['pH']}"

        texts.append(
            ax.text(
                row[x_col],
                row[y_col],
                lbl,
                fontsize=7
            )
        )

    if HAS_ADJUST_TEXT:
        adjust_text(
            texts,
            ax=ax,
            expand_text=(1.2, 1.3),
            expand_points=(1.5, 1.5),
            force_text=0.7,
            force_points=0.5,
            lim=200,
            arrowprops=dict(
                arrowstyle="-",
                lw=0.5,
                alpha=0.6
            )
        )

    ax.set_title(
        title,
        fontsize=20,
        fontweight='bold',
        pad=14
    )

    ax.grid(alpha=0.3, linestyle="--")

    x_pad = (df[x_col].max() - df[x_col].min()) * 0.10
    y_pad = (df[y_col].max() - df[y_col].min()) * 0.10

    ax.set_xlim(
        df[x_col].min() - x_pad,
        df[x_col].max() + x_pad
    )

    ax.set_ylim(
        df[y_col].min() - y_pad,
        df[y_col].max() + y_pad
    )

    method_name = x_col.split("_")[0]

    if method_name == "tSNE":
        method_name = "t-SNE"

    ax.set_xlabel(
        f"{method_name} 1",
        fontsize=18,
        fontweight='bold'
    )

    ax.set_ylabel(
        f"{method_name} 2",
        fontsize=18,
        fontweight='bold'
    )

    ax.tick_params(
        axis='both',
        which='major',
        labelsize=14
    )

    # ---- Panel label moved DOWN so it doesn't hit the suptitle ----
    ax.text(
        -0.08,
        1.02,
        panel_label,
        transform=ax.transAxes,
        fontsize=32,
        fontweight="bold",
        va="bottom",
        ha="right"
    )


# =========================
# MAIN
# =========================
def main():

    print("\n--- Generating 2x2 Master Grid ---\n")

    struct_csv = Path("Structural_4_Metrics_Coordinates.csv")
    omni_csv   = Path("Omni_8_Metrics_Coordinates.csv")

    if not struct_csv.exists() or not omni_csv.exists():
        print("[ERROR] Cannot find required CSV files.")
        return

    df_struct = pd.read_csv(struct_csv)
    df_omni   = pd.read_csv(omni_csv)

    # =========================
    # FIGURE SETUP
    # =========================
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(22, 20)
    )

    # ---------- ROW 1 : STRUCTURAL ----------
    plot_panel(
        axes[0, 0],
        df_struct,
        "tSNE_1",
        "tSNE_2",
        "Structural Features (RMSD, RMSF, Rg, DSSP) — t-SNE",
        "A"
    )

    plot_panel(
        axes[0, 1],
        df_struct,
        "UMAP_1",
        "UMAP_2",
        "Structural Features (RMSD, RMSF, Rg, DSSP) — UMAP",
        "B"
    )

    # ---------- ROW 2 : INTEGRATED ----------
    plot_panel(
        axes[1, 0],
        df_omni,
        "tSNE_1",
        "tSNE_2",
        "Integrated Structural–Ligand Features — t-SNE",
        "C"
    )

    plot_panel(
        axes[1, 1],
        df_omni,
        "UMAP_1",
        "UMAP_2",
        "Integrated Structural–Ligand Features — UMAP",
        "D"
    )

    # =========================
    # MASTER LEGEND
    # =========================
    add_master_legend(fig)

    # =========================
    # MAIN TITLE — pushed higher so it clears panel labels A/B
    # =========================
    fig.suptitle(
        "Global Conformational Landscape of Hexokinase Isoforms Across pH Conditions",
        fontsize=32,
        fontweight='bold',
        y=0.985
    )

    # =========================
    # TIGHT SUBPLOT SPACING
    # =========================
    fig.subplots_adjust(
        top=0.90,      # <-- was 0.93 (gives suptitle breathing room)
        bottom=0.13,
        left=0.055,
        right=0.985,
        wspace=0.14,
        hspace=0.22
    )

    # =========================
    # SAVE OUTPUT
    # =========================
    out_name = "Figure2_Global_Conformational_Landscape.png"

    fig.savefig(
        out_name,
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.15
    )

    fig.savefig(
        out_name.replace(".png", ".tiff"),
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.15
    )

    plt.close(fig)

    print(f"[SUCCESS] Saved figure as:\n{out_name}")


# =========================
# RUN
# =========================
if __name__ == "__main__":
    main()

