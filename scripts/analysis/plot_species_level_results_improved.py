import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# -----------------------------
# Paths
# -----------------------------
PROJECT = "/mnt/c/Users/Aspasia/Desktop/Thesis ΕΚΕΤΑ/species_level_mock_benchmarking"
summary_path = os.path.join(PROJECT, "results/comparative_analysis/species_level_metrics_summary.tsv")
out_dir = os.path.join(PROJECT, "results/comparative_analysis/plots")
os.makedirs(out_dir, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(summary_path, sep="\t")

# -----------------------------
# Cleaner display names
# -----------------------------
tool_display_map = {
    "QIIME2": "QIIME2",
    "DADA2": "DADA2",
    "Kraken2": "Kraken2",
    "Kraken2_direct10": "Kraken2 filtered"
}

dataset_display_map = {
    "dataset_01_miseq_v1v2_mock_dna": "Dataset 1\nSRR3225701",
    "dataset_02_miseq_v1v2_mock_cells_rbb_glycerol": "Dataset 2\nSRR3225702"
}

df["Tool_display"] = df["Tool"].map(tool_display_map).fillna(df["Tool"])
df["Dataset_display"] = df["Dataset"].map(dataset_display_map).fillna(df["Dataset"])
df["Label"] = df["Dataset_display"] + "\n" + df["Tool_display"]

# Tool order
tool_order = {
    "QIIME2": 0,
    "DADA2": 1,
    "Kraken2": 2,
    "Kraken2_direct10": 3
}
df["tool_order"] = df["Tool"].map(tool_order).fillna(99)
df = df.sort_values(["Dataset", "tool_order"]).reset_index(drop=True)

# -----------------------------
# Colors
# -----------------------------
tool_colors = {
    "QIIME2": "#4C78A8",          # blue
    "DADA2": "#F58518",           # orange
    "Kraken2": "#54A24B",         # green
    "Kraken2_direct10": "#B279A2" # purple
}

metric_colors = {
    "Precision": "#4C78A8",
    "Recall": "#F58518",
    "F1_score": "#54A24B"
}

bar_colors = [tool_colors.get(t, "gray") for t in df["Tool"]]

# -----------------------------
# 1. F1-score by tool
# -----------------------------
plt.figure(figsize=(10, 6))
plt.bar(df["Label"], df["F1_score"], color=bar_colors)
plt.title("Species-level F1-score by tool")
plt.ylabel("F1-score")
plt.xlabel("Dataset and tool")
plt.ylim(0, max(df["F1_score"]) + 0.15)
plt.xticks(rotation=35, ha="right")

legend_elements = [
    Patch(facecolor=tool_colors["QIIME2"], label="QIIME2"),
    Patch(facecolor=tool_colors["DADA2"], label="DADA2"),
    Patch(facecolor=tool_colors["Kraken2"], label="Kraken2"),
    Patch(facecolor=tool_colors["Kraken2_direct10"], label="Kraken2 filtered")
]
plt.legend(handles=legend_elements, title="Tool")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "f1_score_by_tool_improved.png"), dpi=300)
plt.close()

# -----------------------------
# 2. Jaccard by tool
# -----------------------------
plt.figure(figsize=(10, 6))
plt.bar(df["Label"], df["Jaccard"], color=bar_colors)
plt.title("Species-level Jaccard similarity by tool")
plt.ylabel("Jaccard similarity")
plt.xlabel("Dataset and tool")
plt.ylim(0, max(df["Jaccard"]) + 0.15)
plt.xticks(rotation=35, ha="right")
plt.legend(handles=legend_elements, title="Tool")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "jaccard_by_tool_improved.png"), dpi=300)
plt.close()

# -----------------------------
# 3. Precision / Recall / F1 grouped plot
# -----------------------------
x = np.arange(len(df))
width = 0.24

plt.figure(figsize=(14, 7))
plt.bar(x - width, df["Precision"], width=width, color=metric_colors["Precision"], label="Precision")
plt.bar(x, df["Recall"], width=width, color=metric_colors["Recall"], label="Recall")
plt.bar(x + width, df["F1_score"], width=width, color=metric_colors["F1_score"], label="F1-score")

plt.title("Species-level precision, recall, and F1-score")
plt.ylabel("Score")
plt.xlabel("Dataset and tool")
plt.ylim(0, 1.0)
plt.xticks(x, df["Label"], rotation=35, ha="right")
plt.legend(title="Metric")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "precision_recall_f1_by_tool_improved.png"), dpi=300)
plt.close()

# -----------------------------
# 4. Predicted vs expected species richness
# -----------------------------
x = np.arange(len(df))

plt.figure(figsize=(12, 7))
plt.bar(x, df["Predicted_species"], color=bar_colors, label="Predicted species")
plt.plot(
    x,
    df["Expected_species"],
    color="black",
    marker="o",
    linestyle="--",
    linewidth=2,
    label="Expected species"
)

plt.title("Predicted versus expected species richness")
plt.ylabel("Number of species")
plt.xlabel("Dataset and tool")
plt.xticks(x, df["Label"], rotation=35, ha="right")

legend_elements_richness = [
    Patch(facecolor=tool_colors["QIIME2"], label="Predicted: QIIME2"),
    Patch(facecolor=tool_colors["DADA2"], label="Predicted: DADA2"),
    Patch(facecolor=tool_colors["Kraken2"], label="Predicted: Kraken2"),
    Patch(facecolor=tool_colors["Kraken2_direct10"], label="Predicted: Kraken2 filtered"),
    Line2D([0], [0], color="black", marker="o", linestyle="--", label="Expected species")
]
plt.legend(handles=legend_elements_richness, loc="upper right")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "predicted_vs_expected_species_improved.png"), dpi=300)
plt.close()

print("Improved plots saved in:")
print(out_dir)
