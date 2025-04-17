## call somatic SV with long-read data using SAVANA. https://github.com/cortes-ciriano-lab/savana1.20
## set --length: Minimum length SV to consider (default=30),
## set --mapq: Minimum MAPQ of reads to consider (default=5),
## set --depth: Minumum number of supporting reads from tumour OR normal to consider variant (default=3)
## set --buffer: Buffer to add when clustering adjacent (non-insertion) potential breakpoints (default=10)
#####

rule call_somatic_cnvs_savana120_cna:
    input:
        tumour_bam="analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumour_bai="analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai="analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai",
        bp_vcf="analysis/svs/savana/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        phased_vcf = lambda wildcards: get_phased_vcf(wildcards.sample_t, wildcards.flowcell, wildcards.mode),
        genome=config['reference']['file'],
    output:
        "analysis/cnvs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}_fitted_purity_ploidy.tsv",
    log:
        "logs/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.cna.log"
    benchmark:
        "benchmarks/savana1.20/{flowcell}.{mode}.{sample_t}.{sample_n}.cna.benchmark.txt"
    envmodules:
        "conda-envs/base"
    params:
        outdir="analysis/cnvs/savana1.20/{flowcell}/{mode}/{sample_t}.{sample_n}",
        blacklist="/mnt/backedup/home/jiaZ/working/data/hg38-blacklist.v2.bed"
    threads: 30
    resources:
        mem = 400,
        walltime = 48
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana1.20
        set -eu

        if [ "$(ls -A {params.outdir})" ]; then   # need to make sure the outdir is empty
            rm -f {params.outdir}/*
        fi
        savana cna \
        --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --phased_vcf {input.phased_vcf} \
        --breakpoints {input.bp_vcf} \
        --blacklist {params.blacklist} \
        --threads {threads} \
        --outdir {params.outdir} &> {log}
        """
