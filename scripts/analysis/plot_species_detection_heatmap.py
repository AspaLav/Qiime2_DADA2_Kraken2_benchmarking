import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# -----------------------------
# Paths
# -----------------------------
PROJECT = os.getcwd()

details_path = os.path.join(
    PROJECT,
    "results",
    "comparative_analysis",
    "species_level_metrics_details.tsv"
)

out_dir = os.path.join(
    PROJECT,
    "results",
    "comparative_analysis",
    "plots"
)

os.makedirs(out_dir, exist_ok=True)

# -----------------------------
# Helper functions
# -----------------------------
def split_species(x):
    if pd.isna(x):
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "nan":
        return []
    return [s.strip() for s in x.split(";") if s.strip()]


# -----------------------------
# Display names
# -----------------------------
tool_map = {
    "Kraken2_direct10": "Kraken2 filtered",
    "Kraken2_filtered": "Kraken2 filtered"
}

dataset_map = {
    "dataset_01_miseq_v1v2_mock_dna": "Dataset 1\nSRR3225701",
    "dataset_02_miseq_v1v2_mock_cells_rbb_glycerol": "Dataset 2\nSRR3225702"
}

tool_order = [
    "QIIME2",
    "DADA2",
    "Kraken2",
    "Kraken2 filtered"
]

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(details_path, sep="\t")
df["Tool"] = df["Tool"].replace(tool_map)

datasets = [d for d in df["Dataset"].unique() if d in dataset_map]

# -----------------------------
# Plot settings
# -----------------------------
cmap = ListedColormap(["#e6e6e6", "#66c2a5"])

fig, axes = plt.subplots(
    1,
    len(datasets),
    figsize=(6.5 * len(datasets), 10),
    constrained_layout=False
)

if len(datasets) == 1:
    axes = [axes]

# -----------------------------
# Create one heatmap per dataset
# -----------------------------
for ax, dataset in zip(axes, datasets):
    sub = df[df["Dataset"] == dataset].copy()

    expected_species = set()
    detected_by_tool = {}

    for _, row in sub.iterrows():
        tool = row["Tool"]

        correct = set(split_species(row.get("Correct_species", "")))
        missed = set(split_species(row.get("Missed_species", "")))

        expected_species.update(correct)
        expected_species.update(missed)

        detected_by_tool[tool] = correct

    species_list = sorted(expected_species)
    present_tools = [tool for tool in tool_order if tool in sub["Tool"].unique()]

    matrix = np.array([
        [
            1 if species in detected_by_tool.get(tool, set()) else 0
            for tool in present_tools
        ]
        for species in species_list
    ])

    ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    # X and Y labels
    ax.set_xticks(np.arange(len(present_tools)))
    ax.set_xticklabels(present_tools, rotation=30, ha="right", fontsize=10)

    ax.set_yticks(np.arange(len(species_list)))
    ax.set_yticklabels(species_list, fontsize=9)

    ax.set_title(dataset_map[dataset], fontsize=13, fontweight="bold", pad=8)
    ax.set_xlabel("Tool", fontsize=11)

    if ax == axes[0]:
        ax.set_ylabel("Expected species", fontsize=11)
    else:
        ax.set_ylabel("")

    # Add cell borders
    ax.set_xticks(np.arange(-0.5, len(present_tools), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(species_list), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Add stronger vertical separator lines between tools
    for x in np.arange(0.5, len(present_tools) - 0.5, 1):
        ax.axvline(x, color="gray", linewidth=1.4, alpha=0.8)

    # Add checkmarks for detected species
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i, j] == 1:
                ax.text(
                    j,
                    i,
                    "✓",
                    ha="center",
                    va="center",
                    fontsize=12,
                    fontweight="bold",
                    color="black"
                )

# -----------------------------
# Legend and title
# -----------------------------
legend_handles = [
    Patch(facecolor="#66c2a5", edgecolor="black", label="Detected"),
    Patch(facecolor="#e6e6e6", edgecolor="black", label="Not detected")
]

fig.suptitle(
    "Detection of expected species by method",
    fontsize=17,
    fontweight="bold",
    y=0.985
)

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.945),
    ncol=2,
    frameon=False,
    fontsize=10
)

# Leave space at the top for title and legend
plt.subplots_adjust(
    top=0.87,
    bottom=0.12,
    left=0.12,
    right=0.98,
    wspace=0.45
)

# -----------------------------
# Save figure
# -----------------------------
out_file = os.path.join(out_dir, "species_detection_heatmap_expected_improved.png")

plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.close()

print("Saved:", out_file)
