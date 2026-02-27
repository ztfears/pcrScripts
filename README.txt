=========== PCR Primer/Probe Workflow =============
The PCR Primer/Probe workflow is used to find primer/probe setups sensitive to a bacterial species and specific against the remainder of its genus

The input for the workflow is a set of complete genome FASTA files for an entire genus, typically found through BV-BRC.org and downloaded from GenBank
These files should also have a prefix added to the filenames to designate the target species/species group and the rest of the genus, where the workflow assumes you labeled using "spec" and "gen" prefixes
	(a useful command for downloading genomes from a CSV of GenBank Accession numbers can be found in usefulCommands.txt)
	(another useful command for adding prefixes to the beginning of your genome files based on the name found in the file text can also be found at usefulCommands.txt)
The output is a series of BoulderIO formatted text files with the primer/probe setups specific for your target.


The current iteration of the workflow comprises of two parts:
- Preprocessing with a Snakemake pipeline
- Generating Primer/Probe setups using Primer3 and some other more hands-on scripts to clean up the output

Included in this repo are the following:
- The Snakemake pipeline with dependent scripts
- Standalone scripts for running each Snakemake step separately for better troubleshooting
- Scripts for formatting input, submitting to, and processing output from Primer3, along with its necessary settings file(s)
- usefulCommands.txt, which contain some helpful Bash commands for downloading genomes and (re)labelling them
- YAML files for setting up the conda envs Snakemake and Primer3 rely on

Instructions on running any script, whether the submission script for Snakemake or others, is written at the top of said script in a comment so you can easily copy, paste, and run it

The Snakemake pipeline assumes you have the following file structure:

targetInputDir/
├─ strains/
│  ├─ genCP000XXXX.fna
│  ├─ genCP000XXXX.fna
│  ├─ specCP000XXXX.fna
...

with all of your prefix-labelled genome files inside the strains file

The Snakemake pipeline consists of the following major steps, with some minor adjustments steps in between:
1. Gene annotation for each genome (Prokka)
2. Pan-genome analysis (Roary)
3. Finding either species-only genes (99% in species; <1% in genus) or genus-wide genes (99% across entire genus) (findTargetGenes.R)
4. Gathering all alleles for said genes from the various genomes (gatherAlleles.sh)
5. Multi-sequence alignment on these multi-allele files (MAFFT)

which ultimately leads to a folder of AFA files ready to be processed and analyzed with Primer3

The Primer3 toolset is used in the following way:
1. Preprocessing the AFA files for each gene and submitting that data to Primer3 with either submit_ampliconMaker.sh (species-only) or submit_ampliconMaker_mask.sh (genus-wide)
2. Post-processing the output from Primer3, sorting by penalty scores to give the top 50 genes and their corresponding primer/probe setups using parse using parsePrimerOutput.R
   a. The output of this post-processing is a set of .csv files (primer/probe setups) and .txt files (sequence consensus for the gene) for your convenience

Ideally in the future we'll have the Primer3 tools integrated into the Snakemake pipeline, but for the time being they're quite simple to implement


A note on the Primer3 melting temperature settings:
At first glance they may seem exceedingly low for standard PCR primer/probe tools but when double-checked with an oligo Tm calculator from IDT, the Tm's come out to ~60C for primers and ~70C for the probes.
The reason why there's such a discrepancy is that Primer3 utilizes only the SantaLucia method for calculating Tm's, while IDT uses the SantaLucia thermodynamic formula along with the Owczarzy salt correction formula. In practice, members of our lab have found IDT's prediction to be fairly accurate so we lowered our target Tm ranges on Primer3 so that when put through IDT, we'd get our ideal Tm's
