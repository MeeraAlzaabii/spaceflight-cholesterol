# Decisions: what changed from the original code, and why

The analysis code comes from the authors' 2024 BMED 430 project. Only changes needed to make it reproducible were made.

- **Paths:** personal `C:/Users/...` paths replaced by paths relative to the project.
- **Data:** downloaded by `scripts/download_data.py` and checked against
  `metadata/checksums.sha256`, instead of by hand.
- **Packages:** installed from `pixi.toml` and `pixi.lock`, instead of
  `BiocManager::install()` in the original script.
- **Settings:** moved from the code to `config/config.yaml`.
- **Sample labels:** taken from the sample names instead of the column order.
- **KEGG:** gene sets downloaded once and saved, instead of queried live on every run.



## Differences from the 2024 results

The 2024 results from the original analysis in R found 42 significant genes, while this reproducible pipeline finds only 4 significant genes. After research, it was found that GeneLab reprocessed the GLDS-241 dataset after the original project: the current files were created on 24 April 2026. Despite that, the main findings match the original analysis, which are the upregulation of Angptl8 in spaceflight, and  Cholesterol metabolism as the enriched pathway.