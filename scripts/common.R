#Shared helpers for the analysis scripts. 
#Run every script from the project folder , so all paths are relative to it.

cfg <- yaml::read_yaml("config/config.yaml")

FIG <- "results/figures"
TAB <- "results/tables"
LOG <- "results/logs"
PROC <- "data/processed"
for (d in c(FIG, TAB, LOG, PROC)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

#Write a tabseparated table without quotes or row names
write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

#Save a ggplot (factoextra plots are ggplots too)
save_gg <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 150, bg = "white")
}

#Placeholder figure, so a step always produces the files Snakemake expects
empty_png <- function(path, message) {
  png(path, width = 900, height = 500, res = 150)
  plot.new()
  text(0.5, 0.5, message)
  invisible(dev.off())
}

#Record package versions for this step
save_session <- function(step) {
  writeLines(capture.output(sessionInfo()), file.path(LOG, paste0(step, "_sessionInfo.txt")))
}