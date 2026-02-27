library(seqinr)
library(dplyr)
library(tidyr)
library(stringr)


args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}

input_file <- args[1]
consensus_threshold <- as.numeric(args[2])

#input_file <- "~/Desktop/pcrProj/eFaecium/mafft/twofer_zwf_alleles.fna.afa"
#input_file <- "twofer_zwf_alleles.fna.afa"
#consensus_threshold <- 1.0

directory <- dirname(input_file)

#strips everything before the first _ and everything after the last underscore
geneID <- str_replace(str_replace(basename(input_file), "^[^_]+_", ""), "^(.+)_[^_]+$", "\\1")

#geneID <- substr(basename(input_file),6,nchar(basename(input_file))-16)

consensus_prop_fn <- function(chars) {
  if (length(chars) == 0) return(NA_real_)
  freqs <- table(chars)
  max(freqs / sum(freqs))
}

pos_consensus_prop <- function(aln_tbl) {
  aln_tbl |>
    mutate(pos = row_number()) |>
    pivot_longer(
      cols = -pos, 
      names_to = "seq_id", 
      values_to = "base") |>
    group_by(pos) |>
    summarise(
      consensus_prop = consensus_prop_fn(base),
    ) |>
    ungroup()
}

msa <- read.alignment(input_file, format = "fasta")
seq_length <- as.numeric(getLength(msa$seq[1]))

spec_consensus_prop <- pos_consensus_prop(as_tibble(t(as.matrix.alignment(msa))))
con <- str_to_upper(consensus(msa, method = "majority"))


nonconserved_regions <- tibble(
  start = numeric(),
  stop = numeric()
)

isHighlightOn <- FALSE
startPos <- numeric()
stopPos <- numeric()


#compare w threshold, mark those below threshold
for(i in 1:seq_length){
  isBelowThreshold <- spec_consensus_prop[i,2] < consensus_threshold
  if(isBelowThreshold & isHighlightOn){
    next
  } else if(isBelowThreshold & !isHighlightOn){
    startPos <- i
    isHighlightOn <- TRUE
    next
  } else if(!isBelowThreshold & isHighlightOn){
    stopPos <- i-1
    nonconserved_regions <- nonconserved_regions |>
      add_row(
        start = startPos,
        stop = stopPos
      )
    isHighlightOn <- FALSE
    next
  } else if(!isBelowThreshold & !isHighlightOn){
    next
  }
}


### can't do brackets, gotta do SEQUENCE_EXCLUDED_REGION=start,length start,length; plus its zero-indexed
###

nonconserved_regions_formatted <- nonconserved_regions |>
  mutate(
    length = stop - start + 1,
    start = start - 1
  ) |> select(start,length) |>
  mutate(
    string = paste0(start,",",length)
  )

sequence_excluded_regions <- paste0("SEQUENCE_EXCLUDED_REGION=",
                                    paste(nonconserved_regions_formatted$string, collapse = " "))
sequence_internal_excluded_regions <- paste0("SEQUENCE_INTERNAL_EXCLUDED_REGION=",
                                    paste(nonconserved_regions_formatted$string, collapse = " "))

#conserved_regions <- conserved_regions |>
#  arrange(desc(stop))

primer3inVec <- con

# for (p in 1:nrow(conserved_regions)) {
#   primer3inVec <- append(primer3inVec, "}", after = conserved_regions[[p,"stop"]])
#   primer3inVec <- append(primer3inVec, "{", after = (conserved_regions[[p,"start"]] - 1))
# }

primer3inString <- paste0(primer3inVec, collapse = "")

### replace all 5xA 5xT 3xc and 3xg with n

primer3inSeq <- str_replace_all(primer3inString, "A{5,}|T{5,}|G{3,}|C{3,}", function(m) str_dup("N", nchar(m)))

primer3inSeq <- str_replace_all(primer3inSeq, "-", "N")

sequence_template <- paste0("SEQUENCE_TEMPLATE=",primer3inSeq)
id <- geneID
id_string <- paste0("SEQUENCE_ID=",id)

primer3input <- paste(id_string, sequence_excluded_regions, sequence_internal_excluded_regions, sequence_template, "=", sep = "\n")

tempString <- c(directory,"/",id,"_input.txt")
outfile <- paste(tempString, collapse="")
cat(primer3input,file=outfile)
