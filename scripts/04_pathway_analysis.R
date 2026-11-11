#Step 4: pathway analysis (GAGE), pathway diagrams (pathview) and KEGG enrichment of the significant genes (clusterProfiler)
#KEGG gene sets are read from the file saved by `pixi run download`
#enrichKEGG() is replaced by enricher() on those same sets (same test)
source("scripts/common.R")
suppressPackageStartupMessages({
  library(gage)
  library(pathview)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggplot2)
})

kegg <- readRDS(cfg$kegg$gene_sets)
kegg_gene_sets <- kegg$kg.sets
res <- read.delim(file.path(TAB, "deseq2_results.tsv"))
significant_genes <- read.delim(file.path(TAB, "significant_genes.tsv"))

to_entrez <- function(ensembl_ids) {
  ids <- sub("\\.\\d+$", "", ensembl_ids)   # drop Ensembl version suffix if any
  suppressMessages(mapIds(org.Mm.eg.db, keys = ids, column = "ENTREZID",
                          keytype = "ENSEMBL", multiVals = "first"))
}

#------------------------- GAGE pathway analysis -------------------------#
log2fc <- res$log2FoldChange
entrez_ids <- to_entrez(res$gene_id)
valid <- !is.na(log2fc) & !is.na(entrez_ids)
log2fc <- log2fc[valid]
names(log2fc) <- entrez_ids[valid]

gage_results <- gage(log2fc, gsets = kegg_gene_sets, ref = NULL, samp = NULL)

sort_pathways <- function(tab) {
  tab <- as.data.frame(tab)
  tab <- tab[!is.na(tab$p.val), , drop = FALSE]
  key <- if ("q.val" %in% colnames(tab)) tab$q.val else tab$p.val
  tab <- tab[order(key), , drop = FALSE]
  data.frame(pathway = rownames(tab), tab, check.names = FALSE)
}
upregulated_sorted <- sort_pathways(gage_results$greater)
downregulated_sorted <- sort_pathways(gage_results$less)
write_tsv(upregulated_sorted, file.path(TAB, "gage_upregulated_pathways.tsv"))
write_tsv(downregulated_sorted, file.path(TAB, "gage_downregulated_pathways.tsv"))

#------------------------- Pathway diagrams (pathview) -------------------#
pv_out <- file.path(FIG, "pathview")
dir.create(pv_out, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$kegg$pathview_dir, recursive = TRUE, showWarnings = FALSE)
kegg_dir <- normalizePath(cfg$kegg$pathview_dir)   # full path, so it still works after setwd()

draw_pathway <- function(sorted, suffix) {
  if (nrow(sorted) == 0) return(invisible(NULL))
  id <- sub(" .*", "", sorted$pathway[1])
  message("Drawing top ", suffix, " pathway: ", sorted$pathway[1])
  old <- setwd(pv_out)          # pathview writes into the working folder
  on.exit(setwd(old))
  tryCatch(
    pathview(gene.data = log2fc, pathway.id = id, species = cfg$kegg$species,
             out.suffix = suffix, kegg.dir = kegg_dir),
    error = function(e) warning("pathview failed for ", id, ": ", conditionMessage(e)))
}
draw_pathway(upregulated_sorted, "upregulated")
draw_pathway(downregulated_sorted, "downregulated")

#------------------------- KEGG enrichment of significant genes ----------#
term2gene <- data.frame(term = sub(" .*", "", rep(names(kegg_gene_sets), lengths(kegg_gene_sets))),
                        gene = unlist(kegg_gene_sets, use.names = FALSE))
term2name <- data.frame(term = sub(" .*", "", names(kegg_gene_sets)),
                        name = sub("^\\S+ ", "", names(kegg_gene_sets)))

entrez_ids_top <- to_entrez(significant_genes$gene_id)
entrez_ids_top <- unique(entrez_ids_top[!is.na(entrez_ids_top)])

kegg_enrichment <- NULL
if (length(entrez_ids_top) > 0) {
  kegg_enrichment <- enricher(gene = entrez_ids_top,
                              TERM2GENE = term2gene, TERM2NAME = term2name,
                              pvalueCutoff = cfg$analysis$enrich_pvalue,
                              qvalueCutoff = cfg$analysis$enrich_qvalue)
}

enrich_tab <- if (is.null(kegg_enrichment)) data.frame() else as.data.frame(kegg_enrichment)
write_tsv(enrich_tab, file.path(TAB, "kegg_enrichment.tsv"))
message("Enriched KEGG pathways: ", nrow(enrich_tab))

bar_file <- file.path(FIG, "kegg_enrichment_barplot.png")
dot_file <- file.path(FIG, "kegg_enrichment_dotplot.png")
n_top <- cfg$analysis$top_pathways
if (nrow(enrich_tab) > 0) {
  save_gg(barplot(kegg_enrichment, showCategory = n_top) +
            ggtitle(paste("Top", n_top, "Enriched KEGG Pathways")), bar_file)
  save_gg(dotplot(kegg_enrichment, showCategory = n_top) +
            ggtitle(paste("Top", n_top, "Enriched KEGG Pathways")), dot_file)
} else {
  empty_png(bar_file, "No enriched KEGG pathways at the chosen cutoffs")
  empty_png(dot_file, "No enriched KEGG pathways at the chosen cutoffs")
}

writeLines(c(paste("KEGG gene sets retrieved:", kegg$retrieved),
             paste("Number of gene sets:", length(kegg_gene_sets))),
           file.path(TAB, "kegg_snapshot_info.txt"))
save_session("04_pathway_analysis")