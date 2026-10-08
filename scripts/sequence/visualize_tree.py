"""Phylogenetic-tree visualization for Hexokinase sequence analysis.

Purpose
-------
Reads the sequence-analysis tree/alignment information expected by the
script and produces a publication-oriented visualization of the
Hexokinase phylogeny.

Inputs
------
The input tree/alignment filenames and labels defined in the script.

Outputs
-------
A Matplotlib-rendered phylogenetic-tree figure saved to the output
location defined in the script.

Run from
--------
Run from the directory containing the sequence-analysis input files
expected by the script.

Usage
-----
``python scripts/sequence/visualize_tree.py``

Dependencies
------------
Python, Biopython, Matplotlib, and the sequence-analysis input files.
"""

from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------
# Load tree
# ---------------------------------------------------

tree = Phylo.read("hexokinases_tree.nwk", "newick")

# Midpoint rooting
tree.root_at_midpoint()

# ---------------------------------------------------
# Clean labels
# ---------------------------------------------------

for clade in tree.find_clades():
    if clade.name:
        clade.name = clade.name.split("|")[-1]

# ---------------------------------------------------
# Colors
# ---------------------------------------------------

colors = {
    "HXK1": "#D98C00",   # orange
    "HXK2": "#4EA3D9",   # blue
    "HXK3": "#1B9E77",   # green
    "HXK4": "#B07AA1",   # purple
    "HKDC1": "#D55E00"   # red-orange
}

# ---------------------------------------------------
# Create figure
# ---------------------------------------------------

fig = plt.figure(figsize=(10, 8), dpi=600)
ax = fig.add_subplot(1, 1, 1)

# ---------------------------------------------------
# Draw tree
# ---------------------------------------------------

Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    show_confidence=True
)

# ---------------------------------------------------
# Clean axes
# ---------------------------------------------------

ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticks([])
ax.set_yticks([])

for spine in ax.spines.values():
    spine.set_visible(False)

# ---------------------------------------------------
# Style labels + bootstrap values
# ---------------------------------------------------

for text in ax.texts:

    label = text.get_text()

    # Bootstrap/confidence values
    try:
        float(label)

        text.set_fontsize(6)
        text.set_color("gray")

        continue

    except:
        pass

    # Default
    color = "black"

    # Assign isoform colors
    for isoform, c in colors.items():
        if isoform in label:
            color = c

    # Style labels
    text.set_color(color)
    text.set_fontsize(11)
    text.set_fontstyle("italic")

# ---------------------------------------------------
# Make branches slightly thicker
# ---------------------------------------------------

for line in ax.get_lines():
    line.set_linewidth(1.5)

# ---------------------------------------------------
# Title
# ---------------------------------------------------

plt.title(
    "Evolutionary Phylogeny of Hexokinase Proteins",
    fontsize=16,
    fontweight="bold"
)

# ---------------------------------------------------
# Legend
# ---------------------------------------------------

legend_elements = [

    Line2D(
        [0], [0],
        color=colors["HXK1"],
        lw=3,
        label="HK1"
    ),

    Line2D(
        [0], [0],
        color=colors["HXK2"],
        lw=3,
        label="HK2"
    ),

    Line2D(
        [0], [0],
        color=colors["HXK3"],
        lw=3,
        label="HK3"
    ),

    Line2D(
        [0], [0],
        color=colors["HXK4"],
        lw=3,
        label="HK4"
    ),

    Line2D(
        [0], [0],
        color=colors["HKDC1"],
        lw=3,
        label="HKDC1"
    )
]

plt.legend(
    handles=legend_elements,
    loc="lower left",
    fontsize=9,
    frameon=False
)

# ---------------------------------------------------
# Layout + Save
# ---------------------------------------------------

plt.tight_layout()

plt.savefig(
    "hexokinase_phylogeny_final.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "hexokinase_phylogeny_final.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

