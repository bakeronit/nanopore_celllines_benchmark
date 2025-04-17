## call somatic SV with long-read data using SAVANA. https://github.com/cortes-ciriano-lab/savana1.20
## set --length: Minimum length SV to consider (default=30),
## set --mapq: Minimum MAPQ of reads to consider (default=5),
## set --depth: Minumum number of supporting reads from tumour OR normal to consider variant (default=3)
## set --buffer: Buffer to add when clustering adjacent (non-insertion) potential breakpoints (default=10)
#####

rule call_somatic_sv_savana120_run:
    input:
        tumour_bam="analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumour_bai="analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai="analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai",
        phased_vcf = lambda wildcards: get_phased_vcf(wildcards.sample_t, wildcards.flowcell, wildcards.mode),
        genome=config['reference']['file'],
    output:
        "analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.bedpe",
        "analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints_read_support.tsv",
        "analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        "analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.inserted_sequences.fa",
        #"analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.sv_breakpoints.vcf",
        #"analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf",
    log:
        "logs/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.run.log"
    benchmark:
        "benchmarks/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.run.benchmark.txt"
    envmodules:
        "conda-envs/base"
    params:
        outdir="analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}"
    threads: 30
    resources:
        mem = 900,
        walltime = 48
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana1.20
        set -eu

        if [ "$(ls -A {params.outdir})" ]; then   # need to make sure the outdir is empty
            rm -f {params.outdir}/*
        fi
        savana run \
        --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --threads {threads} \
        --outdir {params.outdir} &> {log}
        """

rule call_somatic_sv_savana120_classify:
    input:
        "analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
    output:
        vcf="analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.vcf",
        somatic_vcf="analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf",
    log:
        "logs/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.classify.log"
    benchmark:
        "benchmarks/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.classify.benchmark.txt"
    envmodules:
        "conda-envs/base"
    params:
        outdir="analysis/svs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}"
    threads: 8
    resources:
        mem = 10,
        walltime = 24
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana1.20
        set -eu

        savana classify
        --vcf {input} \
        --output {output.vcf} \
        --somatic_output {output.somatic_vcf} \
        --ont \
        --cna_rescue \
        --threads {threads}  &> {log}
        """
