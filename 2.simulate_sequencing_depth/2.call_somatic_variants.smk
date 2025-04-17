#!/usr/bin/env python3

#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
ClairS_sif = config['clairS']['sif']
DP_somatic_sif = "/mnt/backedup/home/jiaZ/working/imgs/deepsomatic/deepsomatic_1.6.0.sif"
#samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
wildcard_constraints:
    sample = "|".join(samples),
    depth = "|".join(["60x","45x","30x","15x"])

include: "helper.smk"
include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/snv_calling/clairS.smk"
include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/snv_calling/deepsomatic.smk"

rule all:
    input:
        #expand("analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz",my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.somatic.vcf.gz",my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"])
        expand("analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz",my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.somatic.vcf.gz",my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"])

use rule call_somatic_snv_clairS as call_somatic_snv_clairS_a with:
    input:
        genome = config['reference']['file'],
        tumor_bam = "analysis/bam/{sample_t}.{depth_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.{depth_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.{depth_n}.bam.bai"
    output:
        "analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz"
    params:
        platform = "ont_r10_dorado_sup_5khz",
        outdir = "analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}",
        clair3_model = "/mnt/backedup/home/jiaZ/working/data/ont_models/clair3_models/r1041_e82_400bps_sup_v430",
        indel_option = "--enable_indel_calling"
    log:
        "logs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}.log"
    benchmark:
        "logs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}.benchmark.txt"

rule call_somatic_snv_deepsomatic_a:
    input:
        genome = config['reference']['file'],
        tumor_bam = "analysis/bam/{sample_t}.{depth_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.{depth_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.{depth_n}.bam.bai"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz"
    params:
        model="/mnt/backedup/home/jiaZ/working/data/ont_models/dpsomatic_model/weights-143-0.987994.ckpt"
    threads: 48
    resources:
        mem=40,
        walltime=100
    log:
        directory("logs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}")
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        singularity exec {DP_somatic_sif} /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=ONT_R104 \
        --ref={input.genome} \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
         --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --num_shards={threads} \
        --logging_dir={log} \
        --customized_model={params.model}
        """

rule filter_somatic:
    input:
        "analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.somatic.vcf.gz"
    envmodules:
        "bcftools/1.19",
        "htslib/1.19.1"
    threads: 1
    resources:
        mem=1,
        walltime=1
    shell:
        """
        bcftools filter -i 'GT="1/1"' {input} |bgzip > {output}
        tabix -p vcf {output}
        """
