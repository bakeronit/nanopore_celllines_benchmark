## call somatic SV with long-read data using SAVANA. https://github.com/cortes-ciriano-lab/savana
## set --length: Minimum length SV to consider (default=30),
## set --mapq: Minimum MAPQ of reads to consider (default=5),
## set --depth: Minumum number of supporting reads from tumour OR normal to consider variant (default=3)
## set --buffer: Buffer to add when clustering adjacent (non-insertion) potential breakpoints (default=10)
#####

rule call_somatic_sv_savana:
    input:
        tumour_bam="analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumour_bai="analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai="analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai",
        genome=config['reference']['file'],
    output:
        "analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.bedpe",
        "analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints_read_support.tsv",
        "analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        "analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.sv_breakpoints.vcf",
        "analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf"
    log:
        "logs/savana/{flowcell}.{mode}.{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/savana/{flowcell}.{mode}.{sample_t}.{sample_n}.benchmark.txt"
    envmodules:
        "conda-envs/base"
    params:
        outdir="analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}"
    threads: 24
    resources:
        mem = 240,
        walltime = 48
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana
        set -eu

        if [ "$(ls -A {params.outdir})" ]; then   # need to make sure the outdir is empty
            rm -f {params.outdir}/*
        fi
        savana --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --threads {threads} \
        --outdir {params.outdir}
        """
