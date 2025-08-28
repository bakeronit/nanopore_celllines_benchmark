"""
Perform the same benchmark but with non-standard chromosome excluded
"""
include: "helper.smk"
include: "../2.simulate_sequencing_depth/helper.smk"

configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/bc.config.yaml"


colo829_samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
hcc1937_samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
rule all:
    input:
        expand("analysis/benchmark/depth/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/indels/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/indels/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/indels/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/depth/indels/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        

def get_vcf(wildcards):
    tool = wildcards.tool
    variant_type = wildcards.variant_type
    sample_t = wildcards.sample_t
    sample_n = wildcards.sample_n
    depth_t = wildcards.depth_t
    depth_n = wildcards.depth_n
    vcf_files = {
        'deepsomatic': {
            'snvs': depth_workdir / 'analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.somatic.vcf.gz',
            'indels': depth_workdir / 'analysis/snvs/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/norm_indel.vcf.gz'
        },
        'clairS': {
            'snvs': depth_workdir / 'analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.vcf.gz',
            'indels': depth_workdir / 'analysis/snvs/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/norm_indel.vcf.gz'
        }
    }
    return vcf_files[tool][variant_type]

def get_goldstandard(wildcards): 
    sample = wildcards.sample
    variant_type = wildcards.variant_type[:-1]
    return config[variant_type]['somatic'][sample]

rule clean_goldstandard:
    input:
        get_goldstandard
    output:
        "analysis/benchmark/goldstandard/{variant_type}/{sample}.clean_chrom.vcf.gz"
    threads: 1
    resources:
        mem = 1,
        walltime = 1
    envmodules:
        "bcftools/1.16",
        "htslib/1.19.1"
    shell:
        """
        if [[ ! -f {input}.tbi ]]; then
            bcftools index -t {input}
        fi
        bcftools view --regions `echo chr{{1..22}} chr{{X,Y}} | tr ' ' ','` {input} | bgzip > {output}
        """

use rule clean_goldstandard as clean_vcf with:
    input:
        get_vcf
    output:
        temp("analysis/benchmark/depth/{variant_type}/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.clean_chrom.vcf.gz")

ClairS_sif = config['clairS']['sif']
rule benchmark_vcf:
    input:
        vcf = rules.clean_vcf.output,
        gs = lambda wildcards: f"analysis/benchmark/goldstandard/{wildcards.variant_type}/{wildcards.sample_t.split('_')[0]}.clean_chrom.vcf.gz",
        genome = config['reference']['file']
    output:
        "analysis/benchmark/depth/{variant_type}/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt"
    params:
        outdir = "analysis/benchmark/depth/{variant_type}/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}"
    threads: 4
    resources:
        mem = 20,
        walltime = 2
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        if [[ {wildcards.variant_type} == 'snvs' ]]; then
            echo -e 'Filter\tType\tPrecision\tRecall\tF1-score\tTP\tFP\tFN' > {output}
            echo -ne 'PASS\t' >> {output}
            singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
            --truth_vcf_fn {input.gs} \
            --input_vcf_fn {input.vcf} \
            --input_filter_tag 'PASS' \
            --threads {threads} \
            --ref_fn {input.genome} \
            --output_dir {params.outdir} --roc_fn {params.outdir}/roc | grep -E 'SNV' >> {output}
        else
            echo -e "Type\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output}
            singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
            --truth_vcf_fn {input.gs} \
            --input_vcf_fn {input.vcf} \
            --input_filter_tag 'PASS' \
            --benchmark_indel \
            --threads {threads} \
            --ref_fn {input.genome} \
            --output_dir {params.outdir} | grep -E 'INDEL|INS|DEL' | awk '{{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7}}' >> {output}
        fi
        """
