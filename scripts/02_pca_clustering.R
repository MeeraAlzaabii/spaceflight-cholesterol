#Step 2: normalisation, PCA, k-means and hierarchical clustering.
source("scripts/common.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(factoextra)
})

input <- readRDS(file.path(PROC, "counts.rds"))
conditions <- as.character(input$samples$condition)

#-------------------------- PCA --------------------------------------#
dge <- DGEList(counts = input$counts, group = conditions)

#Filter out lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

#Normalize using TMM, then log2 counts per million
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

pca_res <- prcomp(t(logCPM), center = TRUE, scale. = TRUE)
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- data.frame(sample = colnames(logCPM),
                     PC1 = pca_res$x[, 1],
                     PC2 = pca_res$x[, 2],
                     Condition = conditions)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA of Normalized RNA-Seq Data") +
  theme_minimal() +
  scale_color_manual(values = c("Ground Control" = "red", "Space Flight" = "blue"))
save_gg(p_pca, file.path(FIG, "pca.png"))

write_tsv(pca_df, file.path(TAB, "pca_coordinates.tsv"))
write_tsv(data.frame(component = paste0("PC", seq_along(percentVar)),
                     percent_variance = percentVar),
          file.path(TAB, "pca_variance.tsv"))
write_tsv(data.frame(genes_before_filter = nrow(input$counts),
                     genes_after_filter = nrow(logCPM)),
          file.path(TAB, "filtering_summary.tsv"))

#-------------------------- K-means Clustering -----------------------#
reduced_data <- pca_df[, c("PC1", "PC2")]
rownames(reduced_data) <- pca_df$sample

set.seed(cfg$analysis$seed)   # the elbow plot also runs k-means
p_elbow <- fviz_nbclust(x = reduced_data, FUNcluster = kmeans, method = "wss") +
  labs(title = "Optimal Number of Clusters (Elbow Method)")
save_gg(p_elbow, file.path(FIG, "kmeans_elbow.png"))

k <- cfg$analysis$kmeans_k
set.seed(cfg$analysis$seed)
kmeans_result <- kmeans(reduced_data, centers = k, nstart = cfg$analysis$kmeans_nstart)

clusters <- data.frame(sample = pca_df$sample, condition = conditions,
                       cluster = kmeans_result$cluster)
write_tsv(clusters, file.path(TAB, "kmeans_clusters.tsv"))
print(table(Condition = conditions, Cluster = kmeans_result$cluster))

p_kmeans <- fviz_cluster(kmeans_result, data = reduced_data,
                         geom = "point", stand = FALSE,
                         ellipse.type = "norm",
                         palette = c("#FC4E07", "#00AFBB"),
                         ggtheme = theme_minimal(),
                         main = "K-Means Clustering (2 Clusters)")
save_gg(p_kmeans, file.path(FIG, "kmeans_clusters.png"))

#-------------------------- Hierarchical Clustering ------------------#
dist_matrix <- dist(reduced_data, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = cfg$analysis$hclust_method)

status_colors <- ifelse(pca_df$Condition == "Ground Control", "#FC4E07", "#00AFBB")

p_dend <- fviz_dend(hclust_result,
                    k = k,
                    label_cols = status_colors[hclust_result$order],  #colours in leaf order
                    k_colors = c("#FC4E07", "#00AFBB"),
                    rect = TRUE,
                    rect_fill = TRUE,
                    main = "Hierarchical Clustering Dendrogram",
                    repel = TRUE)
save_gg(p_dend, file.path(FIG, "hierarchical_dendrogram.png"), width = 9)

message("PC1: ", percentVar[1], "%  PC2: ", percentVar[2], "%")
save_session("02_pca_clustering")