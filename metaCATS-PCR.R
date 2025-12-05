####calculate by base conservation for each strain
### solid primers: all sequence is sensitive/specific to strain
### plated primers: most of sequence is conserved for everyone but 3' tip is sensitive/specific to strain

###solid primer:find dom base among each strain, calculate yj based on that
###find primer that has all positions above a certain yj threshold
library(seqinr)
library(dplyr)
library(tidyr)
library(stringr)

primer_length <- 20

args=(commandArgs(TRUE))
if(length(args)==0){
   print("No arguments supplied.")
}
 
fasta_file <- args[1]
metadata_file <- args[2]


#fasta_file <- "all_nusG_alleles.fna.afa"
#metadata_file <- "all_nusG_locVstr.csv"

solid_yj_threshold <- 0.7
plated_yj_threshold <- 0.9

gene_name <- substr(metadata_file,5, nchar(metadata_file)- 12)
directory <- dirname(fasta_file)

xsqGOF <- function(bases) {
  contingency_table <- table(bases)
  if(length(contingency_table) < 2){
    return(0)
  }
  chi_GOF <- chisq.test(contingency_table)
  chi_GOF$p.value
}

xsqIndep <- function(bases, locus_tags) {
  contingency_table <- table(bases, locus_tags)
  chi_indep <- chisq.test(contingency_table)
  chi_indep$p.value
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

entropy_fn <- function(chars) {
  if (length(chars) == 0) return(NA_real_)
  freqs <- table(chars)
  p <- freqs / sum(freqs)
  -sum(p * log(p, base = 2))
}


# Read metadata
metadata <- as_tibble(read.csv(metadata_file, header = TRUE))

# Read and align FASTA sequences
fasta_data <- read.alignment(fasta_file, format = "fasta")
consensus_all <- consensus(fasta_data, method = "majority")

#species-specific/genus-specific consensus sequences
seqVnam <- tibble(
  locus_tag = fasta_data$nam,
  seq = fasta_data$seq
) |>
  mutate(index = row_number()) |>
  left_join(metadata)


spe_seqs <- seqVnam |>
  filter(
    strain == "spe"
  )
alignedSpe <- list(nb = nrow(spe_seqs), nam = spe_seqs$locus_tag, seq = spe_seqs$seq)
class(alignedSpe) <- "alignment"
consensus_spe <- consensus(alignedSpe, method = "majority")

genseqs <- seqVnam |>
  filter(
    strain == "gen"
  )
alignedGen <- list(nb = nrow(genseqs), nam = genseqs$locus_tag, seq = genseqs$seq)
class(alignedGen) <- "alignment"
consensus_gen <- consensus(alignedGen, method = "majority")


# Convert to character matrix and ensure dimensions match
sequence_length <- as.numeric(getLength(fasta_data$seq[1]))

fasta_tbl <- as_tibble(t(as.matrix.alignment(fasta_data))) |>
  mutate(pos = row_number()) |>
  pivot_longer(
    cols = -pos, 
    names_to = "locus_tag", 
    values_to = "base") |>
  left_join(metadata, join_by(locus_tag))

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
    xsqIndep = xsqIndep(base, strain),
    xsqGOF = xsqGOF(base),
    yjspe = youdensJ(base,strain,"spe",domBaseSpe),
    yjgen = youdensJ(base,strain,"gen",domBaseGen),
    ent = entropy_fn(base)
  ) |>
  ungroup()

### SOLID: looking for sequences that are entirely specific and sensitive for among any strains

solid_spe <- tibble(
  start = numeric(),
  stop = numeric(),
  xsqIndeppvalues = "",
  xsqGOFpvalues = "",
  yjspevalues = "",
  yjgenvalues = ""
)

solid_gen <- tibble(
  start = numeric(),
  stop = numeric(),
  xsqIndeppvalues = "",
  xsqGOFpvalues = "",
  yjspevalues = "",
  yjgenvalues = ""
)

n_possible_primers <- sequence_length -(primer_length-1)
for(p in 1:n_possible_primers){
  sliding_window <- fasta_stats |>
    filter(
      pos >= p & pos < (p+primer_length)
    )
  xsqGOFpvalues <- sliding_window |>
    pull(xsqGOF)
  xsqIndeppvalues <- sliding_window |>
    pull(xsqIndep)
  yjspevalues <- sliding_window |>
    pull(yjspe)
  yjgenvalues <- sliding_window |>
    pull(yjgen)
  
  #sort by yj values for sequence, toss to solid_spe,or gen bins
  
  if(all(yjspevalues > solid_yj_threshold)){
    solid_spe <- solid_spe |>
      add_row(
        start = p,
        stop = p + primer_length - 1,
        xsqIndeppvalues = paste(xsqIndeppvalues, collapse = ", "),
        xsqGOFpvalues = paste(xsqGOFpvalues, collapse = ", "),
        yjspevalues = paste(yjspevalues, collapse = ", "),
        yjgenvalues = paste(yjgenvalues, collapse = ", ")
      )
  }
  if(all(yjgenvalues > solid_yj_threshold)){
    solid_gen <- solid_gen |>
      add_row(
        start = p,
        stop = p + primer_length - 1,
        xsqIndeppvalues = paste(xsqIndeppvalues, collapse = ", "),
        xsqGOFpvalues = paste(xsqGOFpvalues, collapse = ", "),
        yjspevalues = paste(yjspevalues, collapse = ", "),
        yjgenvalues = paste(yjgenvalues, collapse = ", ")
      )
  }
  
}


#### PLATED: Looking for primers that have:
####    super conserved regions across all strains
####    super sens/specif 3' ends, forwards and backwards for species

plated_spe <- tibble(
  start = numeric(),
  stop = numeric(),
  xsqIndeppvalues = "",
  xsqGOFpvalues = "",
  entvalues = "",
  yjspe_first_last = "",
  yjgen_first_last = ""
)

plated_gen <- tibble(
  start = numeric(),
  stop = numeric(),
  xsqIndeppvalues = "",
  xsqGOFpvalues = "",
  entvalues = "",
  yjspe_first_last = "",
  yjgen_first_last = ""
)

n_possible_primers <- sequence_length -(primer_length-1)
for(p in 1:n_possible_primers){
  sliding_window <- fasta_stats |>
    filter(
      pos >= p & pos < (p+primer_length)
    )
  xsqGOFpvalues <- sliding_window |>
    pull(xsqGOF)
  xsqIndeppvalues <- sliding_window |>
    pull(xsqIndep)
  entvalues <- sliding_window |>
    pull(ent)
  yjspevalues <- sliding_window |>
    pull(yjspe)
  yjgenvalues <- sliding_window |>
    pull(yjgen)
  
  #sort by yj values for sequence, toss to solid_spe,2, etc. bins
  if(all(xsqGOFpvalues[2:primer_length] < 0.05) | all(xsqGOFpvalues[1:(primer_length-1)] < 0.05)){
    if(yjspevalues[1] > plated_yj_threshold | yjspevalues[primer_length] > plated_yj_threshold){
      plated_spe <- plated_spe |>
        add_row(
          start = p,
          stop = p + primer_length - 1,
          xsqIndeppvalues = paste(xsqIndeppvalues, collapse = ", "),
          xsqGOFpvalues = paste(xsqGOFpvalues, collapse = ", "),
          entvalues = paste(entvalues, collapse = ", "),
          yjspe_first_last = paste(yjspevalues[1], yjspevalues[primer_length], sep = ", "),
          yjgen_first_last = paste(yjgenvalues[1], yjgenvalues[primer_length], sep = ", ")
        )
    }
    if(yjgenvalues[1] > plated_yj_threshold | yjgenvalues[primer_length] > plated_yj_threshold){
      plated_gen <- plated_gen |>
        add_row(
          start = p,
          stop = p + primer_length - 1,
          xsqIndeppvalues = paste(xsqIndeppvalues, collapse = ", "),
          xsqGOFpvalues = paste(xsqGOFpvalues, collapse = ", "),
          entvalues = paste(entvalues, collapse = ", "),
          yjspe_first_last = paste(yjspevalues[1], yjspevalues[primer_length], sep = ", "),
          yjgen_first_last = paste(yjgenvalues[1], yjgenvalues[primer_length], sep = ", ")
        )
    }
  }
}

###bunch of duplicates between plated bins
###return only unique sequences

plated_spe_unique <- plated_spe |>
  anti_join(plated_gen, join_by(start,stop)) |> rowwise() |>
  mutate(sequence = paste(consensus_spe[start:stop], collapse = "")) |> 
  mutate(hasGaps = str_detect(sequence, "-")) |>
  filter(hasGaps == FALSE) |>
  select(-hasGaps,-xsqIndeppvalues)


tempStringSolid <- c(directory,"/solid_",gene_name,"_primers.csv")
outfileSolid <- paste(tempStringSolid, collapse="")
tempStringPlated <- c(directory,"/plated_",gene_name,"_primers.csv")
outfilePlated <- paste(tempStringPlated, collapse="")

if (nrow(solid_spe) > 0) {
  write.csv(file=outfileSolid, solid_spe, row.names = FALSE)  
}

if (nrow(plated_spe_unique) > 0) {
  write.csv(file=outfilePlated, plated_spe_unique, row.names = FALSE)  
}

