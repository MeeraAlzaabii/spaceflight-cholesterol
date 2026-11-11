# Download the KEGG gene sets for mouse once and save them in data/raw/kegg/.
# The analysis then reads this file instead of asking KEGG on every run.
# Skipped if the file is already there.
suppressPackageStartupMessages(library(gage))
cfg <- yaml::read_yaml("config/config.yaml")$kegg

if (file.exists(cfg$gene_sets)) {
  message("KEGG gene sets already here: ", cfg$gene_sets)
} else {
  dir.create(dirname(cfg$gene_sets), recursive = TRUE, showWarnings = FALSE)
  kegg <- kegg.gsets(species = cfg$species)
  saveRDS(list(kg.sets = kegg$kg.sets,
               species = cfg$species,
               retrieved = as.character(Sys.Date())),
          cfg$gene_sets)
  message("Saved ", length(kegg$kg.sets), " KEGG gene sets to ", cfg$gene_sets)
}