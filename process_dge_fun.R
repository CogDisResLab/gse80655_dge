# Function to generate DGE

library(edgeR)


process_dge <- function(dge, design, num_genes = 10) {
  norm_dge <- calcNormFactors(dge)
  norm_dge <- estimateDisp(norm_dge, design = design)

  fit <- glmQLFit(norm_dge, design)
  qlf <- glmQLFTest(fit, contrast = c(-1, 1))

  table <- qlf$table
  table$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(table),
                         column = "SYMBOL", keytype = "ENSEMBL",
                         multiVals = "first")
  table <- table %>%
    rownames_to_column("Gene_ID") %>%
    dplyr::select(Gene_ID, Symbol, everything())

  top <- topTags(qlf, n = num_genes)

  top_gene_names <- mapIds(org.Hs.eg.db, keys = rownames(top),
                           column = "SYMBOL", keytype = "ENSEMBL",
                           multiVals = "first")

  mat <- dge$counts[rownames(top),]

  logcpm <- cpm(dge, log = TRUE)
  mcpm <- logcpm[rownames(top),]

  output <- list(
    fit = fit,
    qltest = qlf,
    summary = summary(decideTests(qlf)),
    count_matrix = mat,
    logcpm_matrix = mcpm,
    row_values = top_gene_names,
    complete_table = table)

  return(output)
}
