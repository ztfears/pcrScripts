Currently trying to get everything under one pcrWorkflow script
But for the time being we'll just have to run everthing command by command

Before you start anything PLEASE add a 'genus' or 'spec' prefix to corresponding strain files, the workflow works off the assumption that you have done so and it's crucial for processing

Workflow goes like this:
1. prokka
2. roary
3. findTargetGenes
4. gatherAlleles
5. mafft
6. entropy (for analyzing species-only genes)
   OR megaCATS (for whole genus genes)

At the end of findTargetGenes, the files you'll be working with will either have a 'spec_' or 'all_' prefix designating whether the script found species-only genes or is using core genes respectively.

Below are the commands you'll want to run for each step

prokka:
./submit_prokka.sh dirW/strains outputDir

roary:
sbatch roarySlurm.sh prokkaOutputDir path_to_outputDir

findTargetGenes:
./findTargetGenesWrapper.sh path_to_genePresenceAbsence.csv outputDir

gatherAlleles:
./gatherAllelesWrapper.sh findTargetGenesOutputDir prokkaOutputDir outputDir

mafft:
./submit_mafft.sh gatherAllelesOutputDir outputDir

entropy:
./entropyWrapper.sh mafftOutputDir outputDir

megaCATS:




