import os
from pathlib import Path
import pandas as pd

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


FNAs = generate_sample_ids(os.path.join(INPUT_DIR, "strains"))
#print(f"DEBUG: Found {len(FNAs)} files")
#print(f"DEBUG: FNAs = {FNAs}")


rule all:
    input:
        expand(os.path.join(INPUT_DIR, "prokka/{fna}.gff"), fna=FNAs)

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
