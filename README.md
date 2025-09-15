# atgap

A snakemake pipeline for assembling *Arabidopsis thaliana* genomes from PacBio HiFi reads.

## Contents

- [atgap](#atgap)
  - [Contents](#contents)
  - [Overview](#overview)
    - [Pipeline](#pipeline)
    - [Caveat](#caveat)
  - [Installation](#installation)
  - [Data files](#data-files)
    - [Input files](#input-files)
      - [Raw sequence files](#raw-sequence-files)
      - [Reference genome](#reference-genome)
    - [Output files](#output-files)
  - [Usage](#usage)
    - [Set up](#set-up)
    - [Run the pipeline](#run-the-pipeline)
  - [Acknowledgements](#acknowledgements)

## Overview

### Pipeline

The pipeline is based on efforts to optimise long-read assembly for *A. thaliana* based on efforts by Detlef Weigel and Magnus Nordborg's research groups.
For additional background, please see:

> Rabanal, Fernando A., et al. "Pushing the limits of HiFi assemblies reveals centromere diversity between two Arabidopsis thaliana genomes." Nucleic acids research 50.21 (2022): 12309-12327. DOI: [10.1093/nar/gkac1115](https://doi.org/10.1093/nar/gkac1115)

Here are the main steps:

1. Convert raw `.bam` files to `.fastq` (if necessary).
2. Assemble raw reads to contigs with `hifiasm`.
3. Scaffold contigs to a modified version of the TAIR 10 reference assembly, with repetitive regions masked.
4. Extract assembled autosomes.
<!-- 5. Align raw reads to the assembly. -->

### Caveat

This pipeline is optimised for *A. thaliana* in that we assume the genome is small and the donor is inbred (i.e. completely heterozygous).
That allows us to optimise memory usage for the assembler for small genomes and to disable purging duplicates (hifiasm can autodetect when heterozygous haplotypes are assembled as tandem duplications).
Be aware that many **more natural accessions of *A thaliana show residual heterozygosity** than people typically realise.
In that case, expect assembly errors!
If this is the case, you will need to fiddle with the parameters `-f` and `-l` in the rule running `hifiasm`.

## Installation

Clone the repo to your project folder:
```sh
git clone https://github.com/ellisztamas/atgap.git
```

A conda environment is provided to install the necessary dependencies.
```sh
conda create env -f environment.yml
```

## Data files

### Input files

#### Raw sequence files

The pipeline takes raw `.bam` or `.fastq` files as an input.
Place the files you want to process in a directory with no other `.bam` or `.fastq` files present.

Output files with have the same or similar basename as the input file (e.g. `my_accession.bam` will generate an assembly file called `1158_test_scaffolded_autosomes.fasta`).
However, raw sequencing files will typically have ugly default names, so consider renaming them before assembly.

I am informed by someone from Pacific Biosciences that it is not necessary to trim adapters from raw hifi reads.

#### Reference genome

By default the pipeline scaffolded contigs to a version of the TAIR10 genome for accession Col-0, but with repetetive regions masked.
The (gzipped) fasta file will be downloaded when you clone the repo, and unzipped when you run the pipeline.

I tried scaffolding several Swedish genomes to different reference genomes and found that although the outputs varied in size, it wasn't clear what was best.
Your mileage may vary.
If you want to map to a different reference genome, it is probably easiest to change the path to your fasta file in rule `scaffold_contigs` in `atgap.smk`.

### Output files

* `fastq`: Raw fastq files (these will either be converted from `.bam` format or soft-linked from the original directory)
* `hifiasm`: Full output of `hifiasm`. See the [hifiasm documentation](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output) for details.
* `unscaffolded_contigs`: Unscaffolded contigs in `.fasta` format.
* `scaffolded_contigs`:
    * Full output of `ragtag.py scaffold` (see the [RagTag documentation](https://github.com/malonge/RagTag/wiki/scaffold)).
    * This also includes a file ending with `_scaffolded_autosomes.fasta` which is the assembled autosomes excluding unscaffolded contigas, with keys corrected to remove the `_RagTag` suffix. **This is probably the file you are most interested in.**

## Usage

### Set up

The pipeline is designed to submit individual steps to the CLIP cluster via the SLURM scheduler, so you only need to run the command to start the pipeline in a terminal (see below).

However, as the pipeline will take a long time, you should probably run this inside a `tmux` window so that the window stays active, even if your local machine goes to sleep.
([Here](https://www.howtogeek.com/671422/how-to-use-tmux-on-linux-and-why-its-better-than-screen/) is a tutorial on getting started with `tmux`).

### Run the pipeline

Below is an example script to run the pipeline via the SLURM job scheduler.
This assumes that:

* This pipeline is in directory `/path/to/repo/atgap`.
* You have a directory `/path/to/raw_data` containing raw `.bam` or `.fastq` files.
* You want the pipeline to write to `/path/to/output`.

In addition:

* The argument `-j N` tells the pipeline to run N jobs in parallel, assuming there are N input files.
* `--rerun-incomplete` tells the pipeline to start any failed steps from scratch if you encounter and error and need to come back to it.
* `--restart-times 2` tells the pipeline to attempt the assembly or scaffolding steps with twice the memory if these fail initially.

Change these inputs as necessary.

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

snakemake \
    --snakefile 02_library/atgap/atgap.smk \
    --executor slurm \
    --config input_dir=$input_dir \
    --directory $outdir \
    --rerun-incomplete \
    --restart-times 2 \
    -j4
```

## Acknowledgements

Thanks to Fernando Rabanal and colleagues at the MPI Tübingen for optimising
the pipeline for *A. thaliana* and providing the masked `.fasta` file for
scaffolding.