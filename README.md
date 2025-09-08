# atgap
A snakemake pipeline for assembling *Arabidopsis thaliana* genomes from PacBio hifi reads.

Here are the main steps:

1. Convert raw `.bam` files to `.fastq` (if necessary).
2. Assemble raw reads to contigs with `hifiasm`.
3. Scaffold contigs to a modified version of the TAIR 10 reference assembly, with repetitive regions masked.
4. Extract assembled autosomes.
5. Align raw reads to the assembly.

## Installation

Clone the repo to your project folder:
```sh
git clone https://github.com/ellisztamas/atgap.git
```

A conda environment is provided to install the necessary dependencies.
```sh
conda create env -f environment.yml
```

## Usage

Here is an example script to run the pipeline.
Change the paths for the directories to the pipeline, inputs files, and where output files should be written.

```sh
#!/usr/bin/env bash

conda activate atgap

# Path to the directory containing this pipeline
pipeline=/path/to/repo/atgap
# Directory containing raw .bam or .fastq files
input_dir=/path/to/raw_data
# Output directory
output_dir=/path/to/output
mkdir -p $output_dir

# Run the pipeline.
snakemake \
    --snakefile $pipeline/atgap.smk \
    --executor slurm \
    --config input_dir=$input_dir \
    --default-resources \
    --directory $output_dir \
    -j1
```

## Acknowledgements

Thanks to Fernando Rabanal and colleagues at the MPI Tübingen for optimising
the pipeline for *A. thaliana* and providing the masked `.fasta` file for
scaffolding.