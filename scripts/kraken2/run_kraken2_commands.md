# Kraken2 workflow commands

These commands describe the Kraken2 workflow used for species-level taxonomic classification of the two 16S rRNA mock-community datasets.

Kraken2 was used as a k-mer-based read classification method. Unlike QIIME2 and DADA2, Kraken2 classifies reads directly against a reference database and does not first infer ASVs.

## Database

A Kraken2-compatible Greengenes database was used.

Database path used in the analysis:

```bash
data/references/kraken2/greengenes
```

The database was checked using:

```bash
kraken2-inspect --db data/references/kraken2/greengenes | head
```

Large Kraken2 database files are not included in the GitHub repository and must be downloaded or rebuilt separately.

---

## Dataset 1: SRR3225701

Dataset name:

```text
dataset_01_miseq_v1v2_mock_dna
```

Input file:

```text
data/raw/dataset_01_miseq_v1v2_mock_dna/SRR3225701.fastq.gz
```

Kraken2 classification was run using:

```bash
kraken2 \
  --db data/references/kraken2/greengenes \
  --threads 4 \
  --report results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701.report \
  --output results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701.kraken \
  data/raw/dataset_01_miseq_v1v2_mock_dna/SRR3225701.fastq.gz
```

The Kraken2 report was used to extract species-level assignments.

---

## Dataset 2: SRR3225702

Dataset name:

```text
dataset_02_miseq_v1v2_mock_cells_rbb_glycerol
```

Input file:

```text
data/raw/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/SRR3225702.fastq.gz
```

Kraken2 classification was run using:

```bash
kraken2 \
  --db data/references/kraken2/greengenes \
  --threads 4 \
  --report results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702.report \
  --output results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702.kraken \
  data/raw/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/SRR3225702.fastq.gz
```

The Kraken2 report was used to extract species-level assignments.

---

## Species-level extraction

Species-level rows were extracted from the Kraken2 report files. The extracted files were saved as:

```text
results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701_species.tsv
results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702_species.tsv
```

These files contained the species-level Kraken2 predictions used in the unfiltered benchmark.

---

## Direct-read filtering

Because the unfiltered Kraken2 output produced many low-abundance species-level predictions, an additional filtered Kraken2 output was generated.

The filtering rule was:

```text
Direct_reads >= 10
```

The filtered species files were saved as:

```text
results/dataset_01_miseq_v1v2_mock_dna/kraken2/SRR3225701_species_direct10.tsv
results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/kraken2/SRR3225702_species_direct10.tsv
```

The filtered Kraken2 output was included as a separate method in the comparative analysis and labelled as `Kraken2_direct10`.

---

## Notes

The `.kraken` read-level output files and Kraken2 database files are large and are not included in the GitHub repository.

The GitHub repository contains the scripts and small result tables needed to reproduce the comparison, but large raw data and reference database files must be downloaded or generated separately.