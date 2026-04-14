library(tidyverse)
library(seqinr)

topNum <- 20


### Change this variable to be the bacteria directory you're working from (the one with the primers and mafft folders)
### Then just hit Ctrl+Shift+Enter; then after it's done you should have a bunch of consensus_seq.txt and datasheet.csv files

targetDir <- "/home/kbarto10/lactobacillus"




all_best_primers <- tibble()
primerDir <- paste0(targetDir,"/primers")
mafftDir <- paste0(targetDir,"/mafft")
files_list <- list.files(path = primerDir , pattern = "_output.txt", full.names = TRUE)

for (file in files_list){
  data_wide <- read_lines(file) |>
    discard(~ .x == "" | .x == "=") |>
    str_split_fixed("=", 2) |>
    as_tibble(.name_repair = "minimal") |>
    set_names(c("variable", "value")) |>
    mutate(value = parse_guess(value)) |>
    pivot_wider(names_from = variable, values_from = value)
  if(as.numeric(data_wide$PRIMER_PAIR_NUM_RETURNED) > 1){
    row_to_add <- data_wide |> select(SEQUENCE_ID, contains("0")) |> select(contains(c("SEQUENCE", "PENALTY", "TM")))
    all_best_primers <- all_best_primers |> bind_rows(row_to_add)
  }
}

top_amplicons_by_pair_penalty <- all_best_primers |>
  arrange(PRIMER_PAIR_0_PENALTY) |>
  slice_head(n = topNum) |>
  pull(SEQUENCE_ID)

top_paths <- paste0(primerDir,"/" ,top_amplicons_by_pair_penalty, "_output.txt")
for(gene in top_amplicons_by_pair_penalty){
  gene_path <- paste0(primerDir,"/" , gene, "_output.txt")
  datasheet <- read_lines(gene_path) |>
    discard(~ .x == "" | .x == "=") |>
    str_split_fixed("=", 2) |>
    as_tibble(.name_repair = "minimal") |>
    set_names(c("variable", "value")) |>
    mutate(value = parse_guess(value)) |>
    filter(str_detect(variable, "\\d")) |>
    filter(str_detect(variable, "PRODUCT", negate = TRUE)) |>
    filter(str_detect(variable, "TM|SEQUENCE")) |>
    mutate(
      NAME = str_replace(variable, "_[^_]+$", ""),
      col = str_extract(variable, "[^_]+$")
    ) |>
    select(-variable) |>
    pivot_wider(names_from = col, values_from = value)
  
  output_path <- paste0(primerDir,"/" , gene, "_datasheet.csv")
  write_csv(datasheet, output_path)
}

#### if doing genus-wide genes, you'll need to do some filtering so consensus is only from the species and not the genus

for(geneId in top_amplicons_by_pair_penalty){
  afaFile <- list.files(path = mafftDir, pattern = geneId, full.names = TRUE)
  msa <- read.alignment(afaFile, format = "fasta")
  con <- paste0(str_to_upper(consensus(msa, method = "majority")), collapse = "")
  cat(con, file=paste0(primerDir,"/",geneId,"_consensus_seq.txt"))
}


# error_entries <- character()
# 
# for (file in files_list){
#   tryCatch({
#     data_wide <- read_lines(file) |>
#       discard(~ .x == "" | .x == "=") |>
#       str_split_fixed("=", 2) |>
#       as_tibble(.name_repair = "minimal") |>
#       set_names(c("variable", "value")) |>
#       mutate(value = parse_guess(value)) |>
#       pivot_wider(names_from = variable, values_from = value)
#     if(as.numeric(data_wide$PRIMER_PAIR_NUM_RETURNED) > 1){
#       row_to_add <- data_wide |> select(SEQUENCE_ID, contains("0")) |> select(contains(c("SEQUENCE", "PENALTY", "TM")))
#       all_best_primers <- all_best_primers |> bind_rows(row_to_add)
#     }
#   }, error = function(e) {
#     error_entries <<- c(error_entries, file)
#     print(paste("Error with entry:", file))
#   })
# }
# 
# cat(error_entries, file = "errors_entries.txt", sep = "\n")