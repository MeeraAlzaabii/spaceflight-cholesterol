# scripts

Code for this project. To run everything, use `pixi run reproduce` from the project folder .



To download the data (`pixi run download`):
- `download_data.py`: downloads the GLDS-241 counts table and sample table from
  NASA OSDR, makes them read-only, and checks them against `metadata/checksums.sha256`.
- `fetch_kegg_gene_sets.R`: saves the KEGG mouse gene sets once to `data/raw/kegg/`



Analysis steps (run in order by `workflow/Snakefile`):

- `common.R`: shared settings and helpers, loaded by every step.
- `01_prepare_counts.R`: reads the counts table and labels each sample.
- `02_pca_clustering.R`: normalisation, PCA, k-means and hierarchical clustering.
- `03_dge.R`: differential expression with DESeq2, MA and volcano plots, heatmap.
- `04_pathway_analysis.R`: GAGE pathways, pathview diagrams, KEGG enrichment.

All settings (thresholds, seed, file names) are in `config/config.yaml`.