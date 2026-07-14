# QIIME2 workflow commands

These commands describe the QIIME2 workflow used for species-level taxonomic classification of the two 16S rRNA mock-community datasets.

The analysis was performed using QIIME2 2024.10.

## Dataset 1: SRR3225701

Dataset name:

```text
dataset_01_miseq_v1v2_mock_dna
```

The raw single-end FASTQ file was imported into QIIME2 using a manifest file.

```bash
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path data/manifest/dataset_01_manifest.tsv \
  --output-path results/dataset_01_miseq_v1v2_mock_dna/qiime2/demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
```

Read quality was inspected using:

```bash
qiime demux summarize \
  --i-data results/dataset_01_miseq_v1v2_mock_dna/qiime2/demux.qza \
  --o-visualization results/dataset_01_miseq_v1v2_mock_dna/qiime2/demux.qzv
```

Denoising was performed with DADA2 single-end mode:

```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs results/dataset_01_miseq_v1v2_mock_dna/qiime2/demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 300 \
  --o-table results/dataset_01_miseq_v1v2_mock_dna/qiime2/table.qza \
  --o-representative-sequences results/dataset_01_miseq_v1v2_mock_dna/qiime2/rep-seqs.qza \
  --o-denoising-stats results/dataset_01_miseq_v1v2_mock_dna/qiime2/denoising-stats.qza
```

Taxonomy was assigned using a Greengenes2-based classifier:

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier data/references/qiime2/gg2_classifier.qza \
  --i-reads results/dataset_01_miseq_v1v2_mock_dna/qiime2/rep-seqs.qza \
  --o-classification results/dataset_01_miseq_v1v2_mock_dna/qiime2/taxonomy.qza
```

Taxonomy was exported:

```bash
qiime tools export \
  --input-path results/dataset_01_miseq_v1v2_mock_dna/qiime2/taxonomy.qza \
  --output-path results/dataset_01_miseq_v1v2_mock_dna/qiime2/exported_taxonomy_gg2
```

## Dataset 2: SRR3225702

Dataset name:

```text
dataset_02_miseq_v1v2_mock_cells_rbb_glycerol
```

The same workflow was applied to dataset 2:

```bash
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path data/manifest/dataset_02_manifest.tsv \
  --output-path results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
```

```bash
qiime demux summarize \
  --i-data results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/demux.qza \
  --o-visualization results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/demux.qzv
```

```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 300 \
  --o-table results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/table.qza \
  --o-representative-sequences results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/rep-seqs.qza \
  --o-denoising-stats results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/denoising-stats.qza
```

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier data/references/qiime2/gg2_classifier.qza \
  --i-reads results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/rep-seqs.qza \
  --o-classification results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/taxonomy.qza
```

```bash
qiime tools export \
  --input-path results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/taxonomy.qza \
  --output-path results/dataset_02_miseq_v1v2_mock_cells_rbb_glycerol/qiime2/exported_taxonomy_gg2
```














