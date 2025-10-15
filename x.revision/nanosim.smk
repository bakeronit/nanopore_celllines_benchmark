rule all:
    input:
        expand("analysis/fastq/{sample}/training.bam", sample=['COLO829','COLO829_BL'])

rule bam2fastq:
    input:
        "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/bam/R10/sup/{sample}.bam"
    output:
        "analysis/fastq/{sample}.fastq.gz"
    envmodules:
        "samtools/1.17",
        "htslib/1.19.1"
    threads:10
    resources:
        mem=20,
        walltime=10
    shell:
        """
        samtools view -h -s 0.01 {input} | \
        samtools fastq -T* -@{threads} | bgzip > {output}
        """

rule read_analysis:
    input:
        fastq = rules.bam2fastq.output,
        bam = "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/bam/R10/sup/{sample}.bam"
    output:
        multiext("analysis/fastq/{sample}/training",
                 "_aligned_region.pkl",
                 "_aligned_reads.pkl",
                 "_ht_length.pkl",
                 "_besthit.bam",
                 "_match.hist",
                 "_mis.hist",
                 "_del.hist",
                 "_ins.hist",
                 "_first_match.hist",
                 "_error_markov_model",
                 "_ht_ratio.pkl",
                 ".bam",
                 "_match_markov_model",
                 "_model_profile",
                 "_processed.bam",
                 "_unaligned_length.pkl",
                 "_error_rate.tsv",
                 "_strandness_rate")
    params:
        prefix = "analysis/fastq/{sample}/training"
    envmodules:
        "conda-envs/base"
    threads: 30
    resources:
        mem=36,
        walltime=48
    shell:
        """
        set -eu
        conda activate /mnt/lustre/working/lab_nicw/jiaZ/local/micromanba_envs/nanosim
        set +eu
        read_analysis.py genome \
        --read {input.fastq} \
        --g_alnm {input.bam} \
        --output {params.prefix} \
        --fastq \
        --num_threads {threads}
        """

#rule generate_reference:
#    input:
#        genome = "~/working/data/reference/reference.fasta"
#        gs_vcf = ""
#    output:
#        "analysis/reference/reference.alt.fasta"
#    envmodules:
#        "conda-envs/base"

rule simulate_reads:
    input:
        rules.read_analysis.output,
        alt_genome = "analysis/reference/reference.alt.fasta"
    output:
        "analysis/fastq/{sample}/{sample}_simulated_reads.fasta"
    envmodules:
        "conda-envs/base"
    threads: 24
    resources:
        mem=36,
        walltime=48
    shell:
        """
        set -eu
        conda activate /mnt/lustre/working/lab_nicw/jiaZ/local/micromanba_envs/nanosim
        set +eu
        simulate_reads.py genome \
        --ref_g {input.alt_genome} \
        --coverage 30 \
        --dna_type linear \
        --model_prefix {params.prefix} \
        --max_len 1000000 \
        --seed 36 \
        --output {output} \
        --fastq \
        --num_threads {threads}
        """
