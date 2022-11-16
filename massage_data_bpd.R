# MDD Differential Expression

library(stringr)
library(tidyverse)
library(org.Hs.eg.db)

file <- "GSE80655-Count-Data.csv"

metadata <- "GSE80655-Sample-Metadata.csv"

pheno <- read_csv(metadata) %>%
  filter(clinical_diagnosis %in% c("Bipolar Disorder", "Control"))

counts <- read_csv(file) %>%
  dplyr::select(gene_id, all_of(pheno$`Sample-ID`))

count_genes <- mapIds(org.Hs.eg.db, keys = counts$gene_id,
                      column = "SYMBOL", keytype = "ENSEMBL",
                      multiVals = "first") %>%
  enframe() %>%
  filter(!is.na(value)) %>%
  dplyr::rename(gene_id = name,
                symbol = value)

counts_symboled <- counts %>%
  inner_join(count_genes, by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(symbol, everything())

grouped_sections <- pheno %>%
  dplyr::select(`Sample-ID`, brain_region, gender, clinical_diagnosis) %>%
  unique %>%
  group_by(gender, brain_region) %>%
  group_split()

group_names <- grouped_sections %>%
  bind_rows %>%
  group_by(gender, brain_region) %>%
  nest %>%
  unite("name", brain_region, gender) %>%
  dplyr::select(-data) %>%
  pull(name)

names(grouped_sections) <- group_names


grouped_sections %>%
  map2(group_names, ~ write_csv(.x, str_glue("data/{.y}-bpd-metadata.csv"))) %>%
  map(~pull(.x, `Sample-ID`)) %>%
  map(~ dplyr::select(counts_symboled, all_of(c("symbol", .x)))) %>%
  map(~ dplyr::group_by(.x, symbol)) %>%
  map(~ dplyr::summarise(.x, across(where(is.numeric), sum))) %>%
  map2(group_names, ~ write_csv(.x, str_glue("data/{.y}-bpd-counts.csv")))
