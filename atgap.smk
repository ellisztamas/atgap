import os

# FASTA file giving the TAIR10 assembly, with repetitive regions masked as missing.
# tair10_masked = workflow.basedir + "/TAIR10.hard_masked.fa.gz"

# Get input directory from config or command line
input_dir = config.get("input_dir", "input")

# Collect sample names from either .fastq or .bam
fastq_samples, = glob_wildcards(f"{input_dir}/{{sample}}.fastq")
bam_samples,   = glob_wildcards(f"{input_dir}/{{sample}}.bam")

SAMPLES = sorted(set(fastq_samples) | set(bam_samples))


rule all:
    input:
        expand(
            'scaffolded_contigs/{sample}/{sample}_scaffolded_autosomes.fasta',
            # "read_alignment/{sample}_aligned.bam",
            sample=SAMPLES
            )

# rule map_reads_to_assembly:
#     input:
#         genome="scaffolded_contigs/{sample}/{sample}_scaffolded_autosomes.fasta",
#         fastq="fastq/{sample}.fastq",
#     output:
#         "read_alignment/{sample}_aligned.bam",
#     log:
#         "read_alignment/{sample}.log",
#     resources:
#         nodes=1,
#         tasks=1,
#         cpus_per_task=38,
#         mem_mb_per_cpu=4096,
#         runtime=24*60,
#     threads: 38,
#     shell:
#         """
#         bwa index {input.genome}
#         bwa mem -t {threads} {input.genome} {input.fastq} > {output}
#         """

rule extract_autosomes:
    input:
        "scaffolded_contigs/{sample}/ragtag.scaffold.fasta",
    output:
        "scaffolded_contigs/{sample}/{sample}_scaffolded_autosomes.fasta",
    log:
        "scaffolded_contigs/{sample}/{sample}_extract_autosomes.log",
    params:
        chr_names="Chr1_RagTag Chr2_RagTag Chr3_RagTag Chr4_RagTag Chr5_RagTag",
    resources:
        qos='rapid',
        mem_mb=512,
        runtime=5,
    shell:
        """
        samtools faidx {input}
        samtools faidx {input} {params.chr_names} > {output}
        gawk -i inplace '/^>/ {{sub("_RagTag", "", $1)}} 1' {output}
        samtools faidx {output}
        """

rule scaffold_contigs:
    input:
        genome=workflow.basedir + "/TAIR10.hard_masked.fa",
        contigs="unscaffolded_contigs/{sample}_unscaffolded_contigs.fasta",
    output:
        "scaffolded_contigs/{sample}/ragtag.scaffold.fasta",
        "scaffolded_contigs/{sample}/ragtag.scaffold.asm.paf",
    log:
        "scaffolded_contigs/{sample}/{sample}_ragtag.log",
    params:
        prefix = "scaffolded_contigs/{sample}"
    resources:
        qos='rapid',
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
        runtime=60,
    shell:
        """
        # If you have run RagTag before and there are files left over, RagTag seems to
        # use these, even if you use the -w flag to overwrite.
        # To be on the safe side, remove the whole directory.
        # if [ -d {params.prefix} ]; then
        #     echo "Removing existing output directory to avoid RagTag caching things."
        #     rm -r {params.prefix}
        # fi
        
        # Run the scaffolding program.
        ragtag.py scaffold \
            -q 60 \
            -f 10000 \
            -i 0.6 \
            -w \
            --remove-small \
            -o {params.prefix} \
            {input.genome} \
            {input.contigs}
        """

rule gfa_to_fasta:
    input:
        "hifiasm/{sample}.bp.p_ctg.gfa",
    output:
        "unscaffolded_contigs/{sample}_unscaffolded_contigs.fasta",
    resources:
        qos='rapid',
        mem_mb=512,
        runtime=5,
    shell:
        """
        awk '/^S/{{print ">"$2; print $3}}' {input} > {output}
        """

rule hifiasm:
    input:
        "fastq/{sample}.fastq"
    output:
        "hifiasm/{sample}.bp.p_ctg.gfa"
    log: 
        "hifiasm/{sample}.log",
    params:
        outdir = "hifiasm",
        prefix = "hifiasm/{sample}"
    resources:
        qos='medium',
        nodes=1,
        tasks=1,
        cpus_per_task=38,
        mem_mb_per_cpu=lambda wildcards, attempt: 4096 * attempt,
        runtime=24*60,
    threads: 38,
    shell:
        """
        mkdir -p {params.outdir}
        hifiasm \
            -o {params.prefix} \
            -t {threads} \
            -f 0 \
            -l 0 \
            {input}
        """

# If FASTQ exists already, just expose it
rule collect_fastq:
    input:
        f"{input_dir}/{{sample}}.fastq"
    output:
        temp("fastq/{sample}.fastq")
    resources:
        qos='rapid',
        mem_mb=512,
        runtime=5,
    shell:
        "ln -s -f $(realpath {input}) {output}"

# If only BAM exists, convert to FASTQ
rule bam_to_fastq:
    input:
        f"{input_dir}/{{sample}}.bam"
    output:
        temp("fastq/{sample}.fastq")
    resources:
        qos='rapid',
        mem_mb=1048,
        runtime=60,
    shell:
        "samtools fastq {input} > {output}"

rule unzip_ref_genome:
    input:
        workflow.basedir + "/TAIR10.hard_masked.fa.gz",
    output:
        workflow.basedir + "/TAIR10.hard_masked.fa",
    resources:
        qos='rapid',
        mem_mb=1048,
        runtime=60,
    shell:
        """
        gzip -dc {input} > {output}
        """
