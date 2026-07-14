import os
import re
import pandas as pd

PROJECT = "/mnt/c/Users/Aspasia/Desktop/Thesis ΕΚΕΤΑ/species_level_mock_benchmarking"

DATASETS = {
    "dataset_01_miseq_v1v2_mock_dna": {
        "run": "SRR3225701",
        "ground_truth": "data/ground_truth/dataset_01_miseq_v1v2_mock_dna_ground_truth.tsv",
        "qiime2": "results/dataset_01_miseq_v1v2_mock_dna/qiime2/exported_taxonomy_gg2/taxonomy.tsv",
        "dada2": "results/dataset_01_miseq_v1v2_mock_dna/dada2/dada2_taxonomy_gg2_assignTaxonomy.tsv",
        "kraken2": "results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701_species.tsv",
        "kraken2_direct10": "results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701_species_direct10.tsv",
    },
    "dataset_02_miseq_v1v2_mock_cells_rbb_glycerol": {
        "run": "SRR3225702",
        "ground_truth": "data/ground_truth/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol_ground_truth.tsv",
        "qiime2": "results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/exported_taxonomy_gg2/taxonomy.tsv",
        "dada2": "results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/dada2/dada2_taxonomy_gg2_assignTaxonomy.tsv",
        "kraken2": "results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702_species.tsv",
        "kraken2_direct10": "results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702_species_direct10.tsv",
    },
}

# Manual harmonization for obvious taxonomy naming differences.
# We can expand this later after seeing the first output.

SYNONYMS = {
    "phocaeicola vulgatus": "bacteroides vulgatus",
    "phocaeicola_a vulgatus": "bacteroides vulgatus",
    "clostridium_t beijerinckii": "clostridium beijerinckii",
    "enterococcus_h_360604 faecalis": "enterococcus faecalis",
    "porphyromonas_a gingivalis": "porphyromonas gingivalis",
    "bacillus_a cereus": "bacillus cereus",
    "lactobacillus_a gasseri": "lactobacillus gasseri",
    "staphylococcus_a aureus": "staphylococcus aureus",
    "streptococcus_a agalactiae": "streptococcus agalactiae",
    "listeria monocytogenes_b_319916": "listeria monocytogenes",
    "propionibacterium acnes": "cutibacterium acnes",
}


def clean_species_name(name):
    if pd.isna(name):
        return None

    name = str(name).strip()

    # Remove common taxonomy prefixes.
    name = re.sub(r"^[a-z]__", "", name)

    # Convert underscores that are used as separators.
    name = name.replace(";", " ")
    name = re.sub(r"\s+", " ", name)

    # Remove bracket-like extra notation.
    name = name.strip()

    # Lowercase for comparison.
    name = name.lower()

    # Remove species placeholder values.
    if name in {"", "s__", "na", "nan", "none", "unassigned"}:
        return None

    # Remove taxonomy prefixes inside strings if present.
    name = re.sub(r"\b[dpkpcofgs]__", "", name)
    name = re.sub(r"\s+", " ", name).strip()

    # Apply explicit synonym map.
    name = SYNONYMS.get(name, name)

    # If Greengenes-style genus suffix exists, remove suffix where useful.
    # Example: "clostridium_t beijerinckii" -> "clostridium beijerinckii"
    parts = name.split()
    if len(parts) >= 2:
        genus = re.sub(r"_[a-z0-9]+$", "", parts[0])
        species = parts[1]
        name = f"{genus} {species}"

    name = SYNONYMS.get(name, name)
    return name


def extract_species_from_taxon_string(taxon):
    if pd.isna(taxon):
        return None

    taxon = str(taxon)

    # QIIME2/Greengenes format: ...; s__Species name
    pieces = [p.strip() for p in taxon.split(";")]
    species_piece = None
    for p in pieces:
        if p.startswith("s__"):
            species_piece = p.replace("s__", "").strip()

    if species_piece is None or species_piece == "":
        return None

    return clean_species_name(species_piece)


def load_ground_truth(path, run):
    df = pd.read_csv(path, sep="\t")

    df = df[df[run] > 0].copy()
    species = set()

    for _, row in df.iterrows():
        name = f"{row['Genus']} {row['Species']}"
        cleaned = clean_species_name(name)
        if cleaned:
            species.add(cleaned)

    return species


def load_qiime2(path):
    df = pd.read_csv(path, sep="\t")
    taxon_col = "Taxon" if "Taxon" in df.columns else df.columns[1]

    species = set()
    for taxon in df[taxon_col]:
        s = extract_species_from_taxon_string(taxon)
        if s:
            species.add(s)

    return species


def load_dada2(path):
    df = pd.read_csv(path, sep="\t")

    species = set()

    if "Species" in df.columns:
        for _, row in df.iterrows():
            sp = row["Species"]
            if pd.isna(sp):
                continue

            if "Genus" in df.columns and not pd.isna(row.get("Genus")):
                name = f"{row['Genus']} {sp}"
            else:
                name = str(sp)

            cleaned = clean_species_name(name)
            if cleaned:
                species.add(cleaned)
    else:
        raise ValueError(f"No Species column found in {path}")

    return species


def load_kraken2(path):
    df = pd.read_csv(path, sep="\t")
    species = set()

    for name in df["Species"]:
        cleaned = clean_species_name(name)
        if cleaned:
            species.add(cleaned)

    return species


def compute_metrics(expected, predicted):
    tp_set = expected & predicted
    fp_set = predicted - expected
    fn_set = expected - predicted

    tp = len(tp_set)
    fp = len(fp_set)
    fn = len(fn_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    jaccard = tp / len(expected | predicted) if len(expected | predicted) > 0 else 0

    return {
        "Expected_species": len(expected),
        "Predicted_species": len(predicted),
        "True_positives": tp,
        "False_positives": fp,
        "False_negatives": fn,
        "Precision": precision,
        "Recall": recall,
        "F1_score": f1,
        "Jaccard": jaccard,
        "Correct_species": "; ".join(sorted(tp_set)),
        "Missed_species": "; ".join(sorted(fn_set)),
        "Extra_species": "; ".join(sorted(fp_set)),
    }


def main():
    output_dir = os.path.join(PROJECT, "results/comparative_analysis")
    os.makedirs(output_dir, exist_ok=True)

    summary_rows = []
    detail_rows = []

    for dataset_name, paths in DATASETS.items():
        run = paths["run"]

        expected = load_ground_truth(os.path.join(PROJECT, paths["ground_truth"]), run)

        predictions = {
            "QIIME2": load_qiime2(os.path.join(PROJECT, paths["qiime2"])),
            "DADA2": load_dada2(os.path.join(PROJECT, paths["dada2"])),
            "Kraken2": load_kraken2(os.path.join(PROJECT, paths["kraken2"])),
            "Kraken2_direct10": load_kraken2(os.path.join(PROJECT, paths["kraken2_direct10"])),
        }

        for tool, predicted in predictions.items():
            metrics = compute_metrics(expected, predicted)

            row = {
                "Dataset": dataset_name,
                "Run": run,
                "Tool": tool,
            }
            row.update({k: v for k, v in metrics.items() if k not in {
                "Correct_species", "Missed_species", "Extra_species"
            }})
            summary_rows.append(row)

            detail_rows.append({
                "Dataset": dataset_name,
                "Run": run,
                "Tool": tool,
                "Correct_species": metrics["Correct_species"],
                "Missed_species": metrics["Missed_species"],
                "Extra_species": metrics["Extra_species"],
            })

    summary = pd.DataFrame(summary_rows)
    details = pd.DataFrame(detail_rows)

    summary_file = os.path.join(output_dir, "species_level_metrics_summary.tsv")
    details_file = os.path.join(output_dir, "species_level_metrics_details.tsv")

    summary.to_csv(summary_file, sep="\t", index=False)
    details.to_csv(details_file, sep="\t", index=False)

    print("\n=== Species-level benchmark summary ===")
    print(summary.to_string(index=False))

    print("\nSaved:")
    print(summary_file)
    print(details_file)


if __name__ == "__main__":
    main()
