library(seqinr)
library(dplyr)
library(tidyr)
library(stringr)


args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}

yj_threshold <- 0.80

input_file <- args[1]
input_labels <- args[2]
consensus_threshold <- as.numeric(args[3])

#input_file <- "~/Desktop/pcrProj/eCloacae/mafft/all_lsrR_alleles.fna.afa"
#input_labels <- "~/Desktop/pcrProj/eCloacae/labels/all_lsrR_locVstr.csv"
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

domBase<- function(bases){
  base_freq <- table(bases)
  names(base_freq)[which.max(base_freq)]
}

sens <- function(bases, locus_tags, target,dominant_base){
  true_pos <- sum(bases == dominant_base & locus_tags == target)
  false_neg <- sum(bases != dominant_base & locus_tags == target)
  true_pos / (true_pos + false_neg)
}

specif <- function(bases,locus_tags, target,dominant_base){
  false_pos <- sum(bases == dominant_base & locus_tags != target)
  true_neg <- sum(bases != dominant_base & locus_tags != target)
  true_neg / (true_neg + false_pos)
}

youdensJ <- function(bases,locus_tags, target,dominant_base){
  sens(bases,locus_tags, target,dominant_base) + specif(bases,locus_tags, target,dominant_base) - 1
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
labels <- as_tibble(read.csv(input_labels, header = TRUE))


conAll <- consensus(msa, method = "majority")


#species-specific/genus-specific consensus sequences
seqVnam <- tibble(
  locus_tag = msa$nam,
  seq = msa$seq
) |>
  mutate(index = row_number()) |>
  left_join(labels)

spe_seqs <- seqVnam |>
  filter(
    strain == "spe"
  )
alignedSpe <- list(nb = nrow(spe_seqs), nam = spe_seqs$locus_tag, seq = spe_seqs$seq)
class(alignedSpe) <- "alignment"
conSpe <- consensus(alignedSpe, method = "majority")

genseqs <- seqVnam |>
  filter(
    strain == "gen"
  )
alignedGen <- list(nb = nrow(genseqs), nam = genseqs$locus_tag, seq = genseqs$seq)
class(alignedGen) <- "alignment"
conGen <- consensus(alignedGen, method = "majority")

fasta_tbl <- as_tibble(t(as.matrix.alignment(msa))) |>
  mutate(pos = row_number()) |>
  pivot_longer(
    cols = -pos, 
    names_to = "locus_tag", 
    values_to = "base") |>
  left_join(labels, join_by(locus_tag))

fasta_interim_spe <- fasta_tbl |>
  filter(
    strain == "spe"
  ) |>
  group_by(pos)|>
  summarize(
    domBaseSpe = domBase(base)
  )

fasta_interim_gen <- fasta_tbl |>
  filter(
    strain == "gen"
  ) |>
  group_by(pos)|>
  summarize(
    domBaseGen = domBase(base)
  )

fasta_interims <- fasta_tbl |>
  left_join(fasta_interim_spe) |>
  left_join(fasta_interim_gen)


fasta_stats <- fasta_interims |>
  group_by(pos) |>
  summarise(
    consensus_prop = consensus_prop_fn(base),
    yjspe = youdensJ(base,strain,"spe",domBaseSpe),
    yjgen = youdensJ(base,strain,"gen",domBaseGen)
  ) |>
  ungroup() |>
  mutate(
    isExcluded = yjspe < yj_threshold & consensus_prop < consensus_threshold
  )


nonconserved_regions <- tibble(
  start = numeric(),
  stop = numeric()
)

isHighlightOn <- FALSE
startPos <- numeric()
stopPos <- numeric()


# fasta_stats[1,3] <- 0.95
# fasta_stats[2,3] <- 0.95
# fasta_stats[3,3] <- 0.95
# fasta_stats[5,3] <- 0.95
# fasta_stats[6,3] <- 0.95



#compare w threshold, mark those below threshold
for(i in 1:seq_length){
  isBelowThreshold <- fasta_stats[i,5]
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

primer3inVec <- conSpe

species_specific_base_pos <- fasta_stats |>
  filter(
    yjspe > yj_threshold
  ) |>
  pull(pos)

primer3inVec[species_specific_base_pos] <- str_to_upper(primer3inVec[species_specific_base_pos])

primer3inString <- paste0(primer3inVec, collapse = "")

### replace all 5xA 5xT 3xc and 3xg with n

primer3inSeq <- str_replace_all(primer3inString, "A{5,}|T{5,}|G{3,}|C{3,}", function(m) str_dup("N", nchar(m))) |>
  str_replace_all("a{5,}|t{5,}|g{3,}|c{3,}", function(m) str_dup("n", nchar(m))) |>
  str_replace_all("-", "n")

sequence_template <- paste0("SEQUENCE_TEMPLATE=",primer3inSeq)
id <- geneID
id_string <- paste0("SEQUENCE_ID=",id)

primer3input <- paste(id_string, sequence_excluded_regions, sequence_internal_excluded_regions, sequence_template, "=", sep = "\n")

tempString <- c(directory,"/",id,"_input.txt")
outfile <- paste(tempString, collapse="")
cat(primer3input,file=outfile)
