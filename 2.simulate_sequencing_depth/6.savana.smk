### use nanomonsv to call somatic SVs. 
## set --single_bnd --use_racon, so it can identify single breakend sv as well.
## use --qv20 preset parameters for nanopore Q20 chemistry. 
## set --min_indel_size 10 (default 50) to see if it can get smallers sv.maybe not worth it.
## nanomonsv_filter_simple_repeat_svtype: Post filtering of simple repeat and add svtype to the txt results
## nanomonsv_postprocess_sbnd: annotate svs and classify sbnd.
###########

configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
colo829_samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
hcc1937_samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
wildcard_constraints:
    sample = "|".join(colo829_samples + hcc1937_samples),
    depth = "|".join(["60x","45x","30x","15x"])

include: "helper.smk"

rule all:
    input:
        expand("analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.classified.somatic.vcf", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.classified.somatic.vcf", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),


rule call_somatic_sv_savana_a:
    input:
        genome = config['reference']['file'],
        tumour_bam = "analysis/bam/{sample_t}.{depth_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.{depth_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.{depth_n}.bam",
    output:
        "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.sv_breakpoints.bedpe",
        "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.sv_breakpoints_read_support.tsv",
        "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.sv_breakpoints.vcf",
        "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.classified.sv_breakpoints.vcf",
        "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.classified.somatic.vcf",
    params:
        outdir = "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}"
    threads: 24
    resources:
        mem = 240,
        walltime = 24
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana
        set -eu

        if [ "$(ls -A {params.outdir})" ]; then
            rm -f {params.outdir}/*
        fi
        savana --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n} \
        --threads {threads} \
        --outdir {params.outdir}
        """
