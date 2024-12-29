args <- commandArgs(trailingOnly = TRUE)
input <- args[1] # Input directory
output <- args[2] # Output directory
filtered_samples_csv <- file.path(output, "filtered_samples.csv")

# Function to check and validate the sample sheet
sample_sheet <- function(input, output) {
  # List subdirectories (one level deep) containing samples
  sample_dirs <- list.dirs(input, full.names = TRUE, recursive = FALSE)
  
  # Initialize an empty list to store sample information
  samples <- list()
  
  # Loop through each subdirectory to find fastq files
  for (dir in sample_dirs) {
    # Check for files with either .fq.gz or .fastq.gz extensions
    read1 <- list.files(dir, pattern = "_1\\.(fastq|fq)\\.gz$", full.names = TRUE)
    read2 <- list.files(dir, pattern = "_2\\.(fastq|fq)\\.gz$", full.names = TRUE)
    
    # Ensure both read1 and read2 exist
    if (length(read1) > 0 && length(read2) > 0) {
      # Extract sample name from the directory name
      sample_name <- basename(dir)
      
      # Append to the samples list
      samples[[length(samples) + 1]] <- data.frame(
        sample = sample_name,
        read1 = read1,
        read2 = read2,
        stringsAsFactors = FALSE
      )
    } else {
      warning(sprintf("Missing read1 or read2 files in directory: %s. Skipping.", dir))
    }
  }
  
  # Combine all sample information into a single data frame
  if (length(samples) == 0) {
    stop("No valid samples found in the input directories.")
  }
  sample_sheet <- do.call(rbind, samples)
  
  # Print the sample sheet
  print(sample_sheet)
  
  # Validate the sample sheet and save filtered data
  valid_samples <- data.frame()
  for (i in 1:nrow(sample_sheet)) {
    sample <- sample_sheet$sample[i]
    read1 <- sample_sheet$read1[i]
    read2 <- sample_sheet$read2[i]
    
    if (is.na(sample) || sample == "") {
      warning(sprintf("Sample name missing for row %d. Skipping this sample.", i))
      next
    }
    if (!file.exists(read1)) {
      warning(sprintf("Read1 file not found for sample '%s'. Skipping this sample.", sample))
      next
    }
    if (!file.exists(read2)) {
      warning(sprintf("Read2 file not found for sample '%s'. Skipping this sample.", sample))
      next
    }
    # If all checks pass, add to valid samples
    valid_samples <- rbind(valid_samples, sample_sheet[i, ])
  }
  
  # Write valid samples to a CSV file
  if (nrow(valid_samples) > 0) {
    write.csv(valid_samples, filtered_samples_csv, row.names = FALSE)
    cat("Filtered sample sheet saved to:", filtered_samples_csv, "\n")
  } else {
    cat("No valid samples found. No output file generated.\n")
  }
}

# Call the function
sample_sheet(input, output)
