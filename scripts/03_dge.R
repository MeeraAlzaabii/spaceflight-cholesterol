#Step 3: differential gene expression with DESeq2, MA plot, volcano plot, heatmap

source("scripts/common.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

input <- readRDS(file.path(PROC, "counts.rds"))

#DESeq2 needs integer counts: round the RSEM expected counts
rounded_data <- round(input$counts)

metadata <- data.frame(condition = input$samples$condition,
                       row.names = input$samples$sample)

dds <- DESeqDataSetFromMatrix(countData = rounded_data,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

#DESeq2's own summary (alpha = 0.1, its default)
writeLines(capture.output(summary(res)), file.path(TAB, "deseq2_summary.txt"))

#Significant genes: thresholds from config/config.yaml
padj_cut <- cfg$analysis$padj_cutoff
lfc_cut <- cfg$analysis$lfc_cutoff
res_clean <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
significant_genes <- res_clean[res_clean$padj < padj_cut &
                                 abs(res_clean$log2FoldChange) > lfc_cut, ]
message("Significant genes: ", nrow(significant_genes))

res_df <- data.frame(gene_id = rownames(res), as.data.frame(res))
write_tsv(res_df, file.path(TAB, "deseq2_results.tsv"))
sig_df <- data.frame(gene_id = rownames(significant_genes), as.data.frame(significant_genes))
sig_df <- sig_df[order(sig_df$padj), ]
write_tsv(sig_df, file.path(TAB, "significant_genes.tsv"))
saveRDS(dds, file.path(PROC, "dds.rds"))

#MA plot
png(file.path(FIG, "ma_plot.png"), width = 1050, height = 750, res = 150)
plotMA(res, main = "MA Plot", ylim = c(-5, 5))
invisible(dev.off())

#Volcano plot (raw p-values, as in the original)
volcano_df <- data.frame(logFC = res$log2FoldChange, pvalue = res$pvalue)
volcano_df$significant <- volcano_df$pvalue < 0.05 & abs(volcano_df$logFC) > lfc_cut
p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(pvalue), color = significant)) +
  geom_point() + theme_minimal() + labs(title = "Volcano Plot")
save_gg(p_volcano, file.path(FIG, "volcano_plot.png"))

#Heatmap of significant genes (variance-stabilised counts)
heatmap_file <- file.path(FIG, "heatmap_significant_genes.png")
if (nrow(significant_genes) >= 2) {
  vsd <- vst(dds, blind = FALSE)
  norm_counts_significant <- assay(vsd)[rownames(significant_genes), ]
  pheatmap(norm_counts_significant,
           scale = "row",
           show_rownames = FALSE,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           annotation_col = metadata,
           main = "Heatmap of Differentially Expressed Genes",
           filename = heatmap_file, width = 10, height = 7)
} else {
  empty_png(heatmap_file, "Fewer than 2 significant genes: no heatmap")
}

save_session("03_differential_expression")