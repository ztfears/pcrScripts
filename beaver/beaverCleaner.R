library(tidyverse)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}

csvFile <- args[1]

# csvFile <- "~/Desktop/strep.csv"

csvDir <- dirname(csvFile)

# Read the CSV (single column, no header)
df <- read_csv(
  csvFile,
  col_names = "accession"
)

cleaned <- df %>%
# Keep only the first entry in each cell, ignore the rest
	mutate(accession = str_trim(str_remove(accession, ",.*$"))) %>%
# Add ".1" suffix to any entry that doesn't already have it
	mutate(accession = if_else(str_detect(accession, "\\.\\d+$"), accession, paste0(accession, ".1")))

# Save result
write_csv(cleaned, paste0(csvDir,"/temp.csv"), col_names = FALSE)
