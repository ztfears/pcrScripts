#recreate entropy in tidyverse
library(seqinr)
library(dplyr)
library(tidyr)
 
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}
 
input_file <- args[1]
#input_file <- "spec_gph_3_alleles.fna.afa"
directory <- dirname(input_file) 

geneID <- substr(basename(input_file),6,nchar(basename(input_file))-16)

#entropy functions
entropy_fn <- function(chars) {
  if (length(chars) == 0) return(NA_real_)
  freqs <- table(chars)
  p <- freqs / sum(freqs)
  -sum(p * log(p, base = 2))
}

pos_entropy <- function(aln_tbl) {
  aln_tbl |>
    mutate(pos = row_number()) |>
    pivot_longer(
      cols = -pos, 
      names_to = "seq_id", 
      values_to = "base") |>
    group_by(pos) |>
    summarise(
      entropy = entropy_fn(base),
    ) |>
    ungroup()
}

sum_entropy <- function(ent_tbl) {
  ent_tbl |>
    filter(
      pos >= p & pos < p+20
    ) |>
    summarise(total_sum = sum(entropy)) |>
    pull(total_sum)
}


msa <- read.alignment(input_file, format = "fasta")
seq_length <- as.numeric(getLength(msa$seq[1]))

spec_ent <- pos_entropy(as_tibble(t(as.matrix.alignment(msa))))
con <- consensus(msa, method = "majority")



#calculate primer candidate entropy stats
primer_tbl <- tibble(
  start = numeric(), 
  stop = numeric(), 
  ent_sum = numeric(), 
  ent_final3_fwd = numeric(), 
  ent_final_fwd = numeric(),
  ent_final3_rev = numeric(), 
  ent_final_rev = numeric(),
  sequence = "",
  geneID = ""
  )

n_possible_primers <- seq_length-(20-1)

for(p in 1:n_possible_primers){
  primer_ent <- sum_entropy(spec_ent)
  
  forward_final3_ent <- spec_ent |>
    filter(
      pos >= p+17 & pos < p+20
    ) |>
    summarise(total_sum = sum(entropy)) |>
    pull(total_sum)
  
  forward_final_ent <- spec_ent |>
    filter(
      pos == p +19
    ) |>
    pull(entropy)
  
  reverse_final3_ent <- spec_ent |>
    filter(
      pos >= p & pos < p+3
    ) |>
    summarise(total_sum = sum(entropy)) |>
    pull(total_sum)
  
  reverse_final_ent <- spec_ent |>
    filter(
      pos == p
    ) |>
    pull(entropy)
  
  #retrieve sequences and add to corresponding entry
  #add geneId to each row
  primer_tbl <- primer_tbl |>
    add_row(
      start = p,
      stop = p+19,
      ent_sum = primer_ent,
      ent_final3_fwd = forward_final3_ent, 
      ent_final_fwd = forward_final_ent,
      ent_final3_rev = reverse_final3_ent, 
      ent_final_rev = reverse_final_ent,
      sequence = paste0(con[p:(p+19)], collapse = ""),
      geneID = geneID
    )
}

top50_primers <- primer_tbl |>
  slice_min(
    ent_sum,
    n=50
  )

tempString <- c(directory,"/",geneID,"_",50,"-best_",20,"mers",".csv")
outfile <- paste(tempString, collapse="")
write.csv(file=outfile,top50_primers, row.names=FALSE)
