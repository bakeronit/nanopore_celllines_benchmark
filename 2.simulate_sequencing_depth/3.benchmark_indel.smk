

configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"

ClairS_sif = config['clairS']['sif']
#samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]

include: "helper.smk"

rule all:
    input:
        #expand("analysis/benchmark/indels/somatic/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/benchmark/indels/somatic/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"])
        expand("analysis/benchmark/indels/somatic/clairS/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/indels/somatic/deepsomatic/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt", my_dirty_combinator,sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"])


def get_somatic_indel(wildcards):
    if wildcards.tool == "clairS":
        return "analysis/snvs/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/indel.vcf.gz"
    else:
        return "analysis/snvs/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/output.somatic.vcf.gz"
    
rule benchmark_somatic_indel_calling_a:
    input:
        vcf = get_somatic_indel,
        genome = config['reference']['file'],
    output:
        "analysis/benchmark/indels/somatic/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/summary.txt"
    params:
        outdir = "analysis/benchmark/indels/somatic/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}",
        norm_vcf = "analysis/snvs/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/norm_indel.vcf.gz",
        #gs = "/mnt/backedup/home/jiaZ/working/general/goldstandard/vcfs/colo829/merged_normed_isec_indel.goldstandard.vcf.gz"
        gs = "/mnt/backedup/home/jiaZ/working/general/goldstandard/vcfs/hcc1937/merged_normed_isec_indel.goldstandard.vcf.gz"
    threads: 2 
    resources:
        mem = 10,
        walltime = 2
    envmodules:
        "singularity/3.7.1",
        "bcftools/1.19",
        "htslib/1.19.1"
    shell:
        """
        bcftools norm -a  --atom-overlaps . {input.vcf} |bgzip > {params.norm_vcf}
        
        echo -e "Type\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {params.norm_vcf} \
        --input_filter_tag 'PASS' \
        --benchmark_indel \
        --threads {threads} \
        --ref_fn {input.genome} \
        --output_dir {params.outdir} | grep -E 'INDEL|INS|DEL' | awk '{{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7}}' >> {output}
        """

