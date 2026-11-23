#Step 1: read the GeneLab counts table and label each sample

source("scripts/common.R")

#Find the counts table among the downloaded files
files <- list.files(cfg$data$raw_dir, pattern = "\\.csv$", full.names = TRUE)
counts_file <- NULL
for (pattern in cfg$data$counts_patterns) {
  hits <- grep(pattern, basename(files))
  if (length(hits) == 1) { counts_file <- files[hits]; break }
}
if (is.null(counts_file)) {
  stop("Could not pick exactly one counts table in ", cfg$data$raw_dir,
       ". Files: ", paste(basename(files), collapse = ", "),
       ". Adjust counts_patterns in config/config.yaml.")
}
message("Counts table: ", counts_file)

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- counts[!grepl("^ERCC-", rownames(counts)), , drop = FALSE]  # drop spike-ins if present
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

#Label samples from their names (the original script assumed the column order instead, which breaks if GeneLab orders the columns differently)
sample_names <- colnames(counts)
condition <- ifelse(grepl(cfg$data$ground_label, sample_names), "Ground Control",
             ifelse(grepl(cfg$data$flight_label, sample_names), "Space Flight", NA))
if (anyNA(condition)) {
  stop("Cannot tell the condition of: ", paste(sample_names[is.na(condition)], collapse = ", "))
}
samples <- data.frame(sample = sample_names,
                      condition = factor(condition, levels = c("Ground Control", "Space Flight")))
print(table(samples$condition))

saveRDS(list(counts = counts, samples = samples, counts_file = basename(counts_file)),
        file.path(PROC, "counts.rds"))
write_tsv(samples, file.path(PROC, "samples.tsv"))
save_session("01_prepare_counts")