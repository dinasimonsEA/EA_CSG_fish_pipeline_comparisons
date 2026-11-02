# Script to filter the data to remove any poor quality reads

# get wd
path <- getwd()

# Dataset-specific directories
output_dir  <- file.path(path, "Data", "Temp", test_data_name)
rds_dir     <- file.path(output_dir, "R_objects")
filt_dir    <- file.path(output_dir, "04_filterAndTrim")
results_dir <- file.path(path, "Results", test_data_name)

# ensure directories exist
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# read in the lists of cutadapt files for trimming
fnFs.cut <- readRDS(file = file.path(rds_dir, "02_fnFs.cut.rds"))
fnRs.cut <- readRDS(file = file.path(rds_dir, "02_fnRs.cut.rds"))

# set up filepaths for output files
fnFs.filtN <- file.path(filt_dir, basename(fnFs.cut))
fnRs.filtN <- file.path(filt_dir, basename(fnRs.cut))

# run filterAndTrim based on user defined parameters
out <- filterAndTrim(
  fnFs.cut, fnFs.filtN,
  fnRs.cut, fnRs.filtN,
  maxN = 0,
  maxEE = maxEE,
  truncQ = truncQ,
  minLen = minLen,
  truncLen = truncLen,
  multithread = TRUE
)

# save outputs
saveRDS(fnFs.filtN, file = file.path(rds_dir, "04_fnFs.filtN.rds"))
write.table(fnFs.filtN,
            file = file.path(filt_dir, "04_fnFs.filtN.txt"),
            row.names = FALSE, col.names = FALSE)

saveRDS(fnRs.filtN, file = file.path(rds_dir, "04_fnRs.filtN.rds"))
write.table(fnRs.filtN,
            file = file.path(filt_dir, "04_fnRs.filtN.txt"),
            row.names = FALSE, col.names = FALSE)

saveRDS(out, file = file.path(rds_dir, "04_out.rds"))

#added to save input and filtered reads
write.table(
  out,
  file = file.path(results_dir, "04_filtered_reads.tsv"),
  sep = "\t",
  quote = FALSE
)

#test both samples and...
valid_idx <- logical(length(fnFs.filtN))
for (i in seq_along(fnFs.filtN)) {
  cat("Checking:", basename(fnFs.filtN[i]), "\n")
  valid_idx[i] <- tryCatch({
    # Test both files
    plotQualityProfile(fnFs.filtN[i])
    plotQualityProfile(fnRs.filtN[i])
    TRUE
  }, error = function(e) {
    cat(
      "Removing sample:",
      basename(fnFs.filtN[i]),
      "\nReason:",
      e$message,
      "\n\n"
    )
    FALSE
  })
}

# ... keep only valid samples
fnFs.filtN_plot <- fnFs.filtN[valid_idx]
fnRs.filtN_plot <- fnRs.filtN[valid_idx]
cat(
  "Retained", length(fnFs.filtN_plot),
  "of", length(valid_idx),
  "samples\n"
)

# generate quality plots on the quality-trimmed data
pdf(
  file = file.path(results_dir, "04_post_trim_quality_plots.pdf"),
  width = 10,
  height = 10
)

for (i in seq_along(fnFs.filtN_plot)) {
print(plotQualityProfile(fnFs.filtN_plot[i]))
print(plotQualityProfile(fnRs.filtN_plot[i]))
}

dev.off() # closes pdf