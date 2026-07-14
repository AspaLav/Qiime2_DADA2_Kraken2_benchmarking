# Standalone DADA2 workflow for Dataset 1
# Dataset: dataset_01_miseq_v1v2_mock_dna
# Run: SRR3225701

library(dada2)

# Run this script from the project root:
# species_level_mock_benchmarking/

PROJECT <- getwd()

dataset <- "dataset_01_miseq_v1v2_mock_dna"
sample_name <- "SRR3225701"

input_fastq <- file.path(PROJECT, "data", "raw", dataset, paste0(sample_name, ".fastq.gz"))

output_dir <- file.path(PROJECT, "results", dataset, "dada2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

filtered_fastq <- file.path(output_dir, paste0(sample_name, "_filtered.fastq.gz"))

# Reference training set for species-level taxonomy assignment
train_set <- file.path(PROJECT, "data", "references", "dada2", "gg2_2024_09_toSpecies_trainset.fa.gz")

# -----------------------------
# 1. Quality filtering
# -----------------------------

filter_results <- filterAndTrim(
  fwd = input_fastq,
  filt = filtered_fastq,
  truncLen = 300,
  trimLeft = 0,
  maxN = 0,
  maxEE = 2,
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.table(
  filter_results,
  file = file.path(output_dir, "filtering_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# -----------------------------
# 2. Error learning
# -----------------------------

err <- learnErrors(filtered_fastq, multithread = TRUE)

# -----------------------------
# 3. Dereplication
# -----------------------------

derep <- derepFastq(filtered_fastq)
names(derep) <- sample_name

# -----------------------------
# 4. ASV inference
# -----------------------------

dada_result <- dada(derep, err = err, multithread = TRUE)

# -----------------------------
# 5. Sequence table construction
# -----------------------------

seqtab <- makeSequenceTable(dada_result)

# -----------------------------
# 6. Chimera removal
# -----------------------------

seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# -----------------------------
# 7. Save read tracking
# -----------------------------

track <- data.frame(
  sample = sample_name,
  input = filter_results[1, "reads.in"],
  filtered = filter_results[1, "reads.out"],
  denoised = sum(getUniques(dada_result)),
  nonchim = sum(seqtab_nochim)
)

write.table(
  track,
  file = file.path(output_dir, "read_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 8. Export ASV table
# -----------------------------

asv_sequences <- colnames(seqtab_nochim)
asv_ids <- paste0("ASV", seq_along(asv_sequences))

asv_table <- as.data.frame(t(seqtab_nochim))
asv_table$Feature.ID <- asv_ids
asv_table$Sequence <- asv_sequences

asv_table <- asv_table[, c("Feature.ID", "Sequence", sample_name)]

write.table(
  asv_table,
  file = file.path(output_dir, "asv_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 9. Export ASV sequences as FASTA
# -----------------------------

fasta_file <- file.path(output_dir, "asv_sequences.fasta")

fasta_lines <- as.vector(rbind(
  paste0(">", asv_ids),
  asv_sequences
))

writeLines(fasta_lines, fasta_file)

# -----------------------------
# 10. Save RDS object
# -----------------------------

saveRDS(seqtab_nochim, file = file.path(output_dir, "seqtab_nochim.rds"))

# -----------------------------
# 11. Taxonomy assignment
# -----------------------------

taxa <- assignTaxonomy(
  seqs = asv_sequences,
  refFasta = train_set,
  multithread = TRUE
)

taxa_df <- as.data.frame(taxa)
taxa_df$Feature.ID <- asv_ids
taxa_df$Sequence <- asv_sequences

taxa_df <- taxa_df[, c("Feature.ID", "Sequence", setdiff(colnames(taxa_df), c("Feature.ID", "Sequence")))]

write.table(
  taxa_df,
  file = file.path(output_dir, "dada2_taxonomy_gg2_assignTaxonomy.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("DADA2 workflow completed for Dataset 1:", dataset, "\n")