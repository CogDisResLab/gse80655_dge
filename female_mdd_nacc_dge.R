# Female nAcc DGE

library(tidyverse)
library(edgeR)
source("process_dge_fun.R")

counts_file <- "data/nAcc_female-mdd-counts.csv"
metadata_file <- "data/nAcc_female-mdd-metadata.csv"

counts <- read_csv(counts_file) %>%
  column_to_rownames("symbol")
metadata <- read_csv(metadata_file) %>%
  select(`Sample-ID`, clinical_diagnosis) %>%
  mutate(diag = if_else(clinical_diagnosis == "Control", "CTL", "MDD"))

groups <- metadata$diag
design_matrix <- model.matrix(~ 0 + groups)

dge <- DGEList(counts = counts, group = groups)

keep <- filterByExpr(dge)
dge_filtered <- dge[keep, , keep.lib.sizes=FALSE]

out_filtered <- process_dge(dge_filtered, design_matrix, num_genes = 1000)

out_filtered$complete_table %>%
  rownames_to_column("gene_name") %>%
  write_csv("results/female_mdd_nacc_dge.csv")

save(out_filtered, file = "results/female_mdd_nacc_dge.RData")
