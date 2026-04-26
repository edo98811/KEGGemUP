suppressPackageStartupMessages({
  library("macrophage")
  library("org.Hs.eg.db")
  library("SummarizedExperiment")
  library("AnnotationDbi")
  library("clusterProfiler")
  library("limma")
  library("igraph")
  library("edgeR")
})

# Load the macrophage dataset ---------------------------------------------------
data(gse, package = "macrophage")
rownames(gse) <- gsub("\\..*", "", rownames(gse))

# limma analysis ---------------------------------------------------------------
condition <- factor(colData(gse)[, "condition_name"])
cell_line <- factor(colData(gse)[, "line"])
design <- model.matrix(~0 + condition + cell_line)

contrast_matrix <- makeContrasts(
  IFNg_vs_naive = conditionIFNg - conditionnaive,
  levels = design
)

dge <- DGEList(assay(gse))
dge <- calcNormFactors(dge)

cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop, ]

voom_mat <- voom(dge, design, plot = FALSE)
fit <- lmFit(voom_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  # Empirical Bayes moderation

# Gene annotation ---------------------------------------------------------------

anns <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(gse),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL",
  multiVals = "first"
)

# limma results ---------------------------------------------------------------
res_macrophage_IFNg_vs_naive_limma <- topTable(
  fit2,
  coef = "IFNg_vs_naive",
  adjust = "fdr",
  number = Inf,
  confint = TRUE
)

res_macrophage_IFNg_vs_naive_limma$ENTREZID <- anns$ENTREZID[match(rownames(res_macrophage_IFNg_vs_naive_limma), anns$ENSEMBL)]
res_macrophage_IFNg_vs_naive_limma$gene_name <- anns$SYMBOL[match(rownames(res_macrophage_IFNg_vs_naive_limma), anns$ENSEMBL)]

# Enrichment analysis ----------------------------------------------------------
de_entrez_IFNg_vs_naive_genes <- anns$ENTREZID[
  (!is.na(res_macrophage_IFNg_vs_naive_limma$adj.P.Val)) &
    (res_macrophage_IFNg_vs_naive_limma$adj.P.Valj <= 0.05)
]

saveRDS(res_macrophage_IFNg_vs_naive_limma, "inst/extdata/limma_res_macrophage.RDS", compress = "xz")
