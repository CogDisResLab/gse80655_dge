# Function to generate DGE

library(edgeR)


process_dge <- function(dge, design, num_genes = 10) {
  norm_dge <- calcNormFactors(dge)
  norm_dge <- estimateDisp(norm_dge, design = design)

  fit <- glmQLFit(norm_dge, design)
  qlf <- glmQLFTest(fit, contrast = c(-1, 1))

  table <- qlf$table
  top <- topTags(qlf, n = num_genes)

  mat <- dge$counts[rownames(top),]

  logcpm <- cpm(dge, log = TRUE)
  mcpm <- logcpm[rownames(top),]

  output <- list(
    fit = fit,
    qltest = qlf,
    summary = summary(decideTests(qlf)),
    count_matrix = mat,
    logcpm_matrix = mcpm,
    row_values = rownames(top),
    complete_table = table)

  return(output)
}
