library(dplyr)
library(tidyr)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}

outputDir <- args[2]

gpa <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE) |>
#gpa <- read.csv("eFaeciumOut2.csv") |>
  select(Gene,ends_with(".fna"))

#calculate number of species strains and number of total genus strains
nStrains <- (ncol(gpa) - 1)
nSpecies <- ncol(gpa |> select(starts_with("spec")))

gpa <- gpa |>
  pivot_longer(
    cols = ends_with("fna"),
    names_to = "strain",
    values_to = "locus_tag"
  ) |>
  filter(locus_tag != "") |>
  mutate(
    strain = substr(strain,1,3)
  )

#gather all species-only genes
specOnlyGenes <- gpa |>
  summarise(
    n=n(), #(ncol(strains)-1)
    .by = c(Gene,strain)
  ) |> 
  ungroup() |>
  complete(Gene,strain, fill = list(n = 0)) |>
  pivot_wider(
    names_from = "strain",
    values_from = "n"
  ) |>
  mutate(
    spe = spe/nSpecies
  ) |>
  filter(
    gen == 0 & spe > 0.99
  )

if(nrow(specOnlyGenes) > 0){
  genesToSave <- specOnlyGenes[,1]
  prefix <- "/spec_"
  
} else {
  #gather all core genes for genus
  coreGenes <- gpa |>
    summarise(
      prop = n()/nStrains,
      .by = c(Gene)
    ) |>
    filter(
      prop > 0.99
    )
  
  genesToSave <- coreGenes[,1]
  prefix <- "/all_"
}

#geneID <- "rplK"

for (geneID in t(genesToSave)){
  #pull all locus tags for gene, put in single column df
  locus_tags <- gpa |>
    filter(
      Gene == geneID
    )
  locus_tags <- locus_tags[,3]
  colnames(locus_tags) <- NULL

  locusVstrain <- gpa |>
    filter(Gene == geneID) |>
    select(strain,locus_tag)

  write.csv(locus_tags, paste0(outputDir,prefix,geneID,".csv"), row.names = FALSE, quote = FALSE)
  write.csv(locusVstrain, paste0(outputDir,prefix,geneID,"_locVstr.csv"), row.names = FALSE, quote = FALSE)
}
