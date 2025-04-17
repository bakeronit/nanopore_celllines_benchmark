### ================================================================
### perform methyaltion frequency calculation for R9 data.
### 1. get fastq files, and get fast5 at the same time
### 2. index the fastq files with fast5
### 3. infer methylation frequecy with bam file and indexed fastq
### 4. calculate methylation frequency and merge from multiple runs
### ================================================================

workdir: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/work"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/snakemake-workflow/config/config.yaml"
F5C="/mnt/backedup/home/jiaZ/working/local/f5c/v1.3/f5c_x86_64_linux"
F5C_GPU="/mnt/backedup/home/jiaZ/working/local/f5c/v1.3/f5c_x86_64_linux_cuda"

rule all:
    input:
        expand("analysis/mod/f5c/R10/COLO829/{run}.meth.tsv", run=['PAK71894','PAK71995']),
        expand("analysis/mod/f5c/R10/COLO829_BL/{run}.meth.tsv", run=['PAK72039','PAK72550'])


rule bam_to_fastq:
    input:
        "analysis/ubam/{flowcell}/sup/{sample}/{run}.ubam"
    output:
        temp("analysis/fastq/{flowcell}/{sample}/{run}.fastq.gz")
    envmodules:
        "samtools/1.17",
        "htslib/1.17"
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    shell:
        """
        samtools fastq -@{threads} {input} | bgzip > {output}
        """

rule index_fastq_with_fast5:
    input:
        fastq = "analysis/fastq/R9/{sample}/{run}.fastq.gz",
        fast5_dir = ""
    output:
        temp("analysis/fastq/R9/{sample}/{run}.fastq.gz.index"),
        temp("analysis/fastq/R9/{sample}/{run}.fastq.gz.index.fai"),
        temp("analysis/fastq/R9/{sample}/{run}.fastq.gz.index.gzi"),
        temp("analysis/fastq/R9/{sample}/{run}.fastq.gz.index.readdb")
    threads: 12
    resources:
        mem = 20,
        walltime = 12
    shell:
        """
        {F5C} index -t {threads} --iop 2 -d {fast5_dir} {fastq}
        """

rule index_fastq_with_blow5:
    input:
        fastq = "analysis/fastq/R10/{sample}/{run}.fastq.gz",
        blow5 = "raw_pod5/temp_blow5/R10/{sample}/{run}.blow5"
    output:
        "raw_pod5/temp_blow5/R10/{sample}/{run}.blow5.idx",
        temp("analysis/fastq/R10/{sample}/{run}.fastq.gz.index"),
        temp("analysis/fastq/R10/{sample}/{run}.fastq.gz.index.fai"),
        temp("analysis/fastq/R10/{sample}/{run}.fastq.gz.index.gzi"),
    threads: 12
    resources:
        mem = 24,
        walltime = 12
    shell:
        """
        {F5C} index -t {threads} --slow5 {input.blow5} {input.fastq}
        """

rule call_methylation:
    input:
        "analysis/bam/{flowcell}/sup/{sample}/{run}.bam.csi",
        "analysis/fastq/{flowcell}/{sample}/{run}.fastq.gz.index",
        "analysis/fastq/{flowcell}/{sample}/{run}.fastq.gz.index.fai",
        "analysis/fastq/{flowcell}/{sample}/{run}.fastq.gz.index.gzi",
        fastq = "analysis/fastq/{flowcell}/{sample}/{run}.fastq.gz",
        ref = config['reference']['file'],
        bam = "analysis/bam/{flowcell}/sup/{sample}/{run}.bam",
        blow5 = "raw_pod5/temp_blow5/{flowcell}/{sample}/{run}.blow5",
    output:
        "analysis/mod/f5c/{flowcell}/{sample}/{run}.meth.tsv"
    threads: 32
    resources:
        mem = 120,
        walltime = 72,
        gpu = 1
    params:
        pore = lambda w: "r9" if w.flowcell=="R9" else "r10"
    shell:
        """
        {F5C_GPU} call-methylation -t {threads} -r {input.fastq} -b {input.bam} -g {input.ref} --slow5 {input.blow5} \
        --secondary=no --pore {params.pore} -x hpc > {output}
        """
