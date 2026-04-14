import os
from pathlib import Path
import pandas as pd
import glob

INPUT_DIR = config.get("input_dir")
SCRIPTS_DIR = config.get("scripts_dir")

def generate_sample_ids(input_dir):
    input_path = Path(input_dir)
    
    files = [f for f in input_path.iterdir() if f.is_file()]
    
    if not files:
        print(f"Warning: No files found in {input_dir}")
        df = pd.DataFrame(columns=["filename"])
    else:
        df = pd.DataFrame({"filename": [f.name for f in files]})
    
    return df["filename"].tolist()

def get_target_genes(wildcards):
    checkpoint_output = checkpoints.findTargetGenes.get().output.outdir
    
    csv_files = glob.glob(os.path.join(checkpoint_output, "*.csv"))
    
    # Extract item names (file stems without .csv extension)
    genes = [Path(f).stem for f in csv_files]    
    return expand(os.path.join(INPUT_DIR, "mafft/{gene}_alleles.fna.afa"), gene=genes)


FNAs = generate_sample_ids(os.path.join(INPUT_DIR, "strains"))
#print(f"DEBUG: Found {len(FNAs)} files")
#print(f"DEBUG: FNAs = {FNAs}")

rule all:
    input:
        os.path.join(INPUT_DIR, "success.txt")

rule success_flag:
    input:
        get_target_genes
    output:
        os.path.join(INPUT_DIR, "success.txt")
    resources:
        runtime=2
    threads: 1
    shell:
        """
        echo "Success.txt" > {output}
        """

rule mafft:
    input:
        os.path.join(INPUT_DIR, "alleles/{gene}_alleles.fna")
    output:
        os.path.join(INPUT_DIR, "mafft/{gene}_alleles.fna.afa")
    conda:
        "/nobackup/archive/grp/grp_Assembly/shared-conda-envs/mafft"
    resources:
        runtime=480,
        mem_mb_per_cpu=2000
    threads: 8
    params:
        input_dir=INPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    shell:
        """
        {params.scripts_dir}/mafft.sh {input} {params.input_dir}/mafft
        """

rule gatherAlleles:
    input:
        csv=os.path.join(INPUT_DIR, "targetGenes/{gene}.csv"),
        tempFile=os.path.join(INPUT_DIR, "alleles/all.ffn")
    output:
        os.path.join(INPUT_DIR, "alleles/{gene}_alleles.fna")
    conda:
        "/grphome/grp_Assembly/seqtk"
    resources:
        runtime=30
    threads: 2
    params:
        input_dir=INPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    shell:
        """
        {params.scripts_dir}/gatherAlleles.sh {input.csv} {input.tempFile} {params.input_dir}/alleles
        """

checkpoint findTargetGenes:
    input:
        os.path.join(INPUT_DIR, "roary/gene_presence_absence.csv")
    output:
        outdir=directory(os.path.join(INPUT_DIR, "targetGenes"))
    conda:
        "/grphome/grp_Assembly/tidy"
    resources:
        runtime=60
    threads: 4
    params:
        input_dir=INPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    shell:
        """
        mkdir -p {params.input_dir}/targetGenes
        Rscript --vanilla {params.scripts_dir}/findTargetGenes.R {input} {params.input_dir}/targetGenes
        """

rule combineFFNs:
    input:
        expand(os.path.join(INPUT_DIR, "prokka/{fna}.ffn"), fna=FNAs)
    output:
        temp(os.path.join(INPUT_DIR, "alleles/all.ffn"))
    resources:
        runtime=15
    threads: 1
    params:
        input_dir=INPUT_DIR
    shell:
        """
        cat {params.input_dir}/prokka/*.ffn > {output}
        """

rule roary:
    input:
        expand(os.path.join(INPUT_DIR, "prokka/{fna}.gff"), fna=FNAs)
    output:
        os.path.join(INPUT_DIR, "roary/gene_presence_absence.csv")
    conda:
        "/nobackup/archive/grp/grp_Assembly/shared-conda-envs/roary"
    resources:
        runtime=600,
        mem_mb_per_cpu=3200
    threads: 20
    params:
        input_dir=INPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    shell:
        """
        {params.scripts_dir}/roary.sh {params.input_dir}/prokka {params.input_dir}/roary
        mv {params.input_dir}/roary*/* {params.input_dir}/roary
        """

rule prokka:
    input:
        os.path.join(INPUT_DIR, "strains/{fna}")
    output:
        os.path.join(INPUT_DIR, "prokka/{fna}.gff"),
        os.path.join(INPUT_DIR, "prokka/{fna}.ffn")
    conda:
        "/nobackup/archive/grp/grp_Assembly/shared-conda-pkgs/envs/prokka"
    resources:
        runtime=15
    threads: 4
    params:
        input_dir=INPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    shell:
        """
        {params.scripts_dir}/prokka.sh {input} {params.input_dir}/prokka
        """
