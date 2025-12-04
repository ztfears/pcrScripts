library(dplyr)
library(tidyr)

#returns genes found in >99% of the species and 1-2% of the genus instead of core genes when no species-only genes are found


args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}
 
genusThreshold <- 0.15
speciesThreshold <- 0.85

gpa <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE) |>
#gpa <- read.csv("sAgalactiaeOut.csv") |>
  select(Gene,ends_with(".fna"))

outputDir <- args[2]
# outputDir <- "C:/Users/zachf/RProjects/roaryOut"

#calculate number of species strains and number of total genus strains
nStrains <- (ncol(gpa) - 1)
nSpecies <- ncol(gpa |> select(starts_with("spec")))
nGenus <- ncol(gpa |> select(starts_with("gen")))

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
    spe = spe/nSpecies,
    gen = gen/nGenus
  ) |>
  filter(
    gen < genusThreshold & spe > speciesThreshold
  )

genesToSave <- specOnlyGenes[,1]
prefix <- "/twofer_"


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
