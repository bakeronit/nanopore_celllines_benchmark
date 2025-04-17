configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
colo829_samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
hcc1937_samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]

wildcard_constraints:
    sample = "|".join(colo829_samples + hcc1937_samples),
    sample_t = "|".join(colo829_samples[:-1] + hcc1937_samples[:-1]),
    sample_n = "|".join(['COLO829_BL','HCC1937_BL']),
    depth = "|".join(["60x","45x","30x","15x"])

include: "helper.smk"

rule all:
    input:
        expand("analysis/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/somatic_SVs/severus_somatic.vcf", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/somatic_SVs/severus_somatic.vcf", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),

Clair3_sif = config['clair3']['sif']

rule call_somatic_sv_severus_a:
    input:
        genome = config['reference']['file'],
        hp_tagged_tumour_bam = "analysis/bam/{sample_t}.{depth_t}.haplotagged.bam",
        hp_tagged_tumour_bai = "analysis/bam/{sample_t}.{depth_t}.haplotagged.bam.bai",
        hp_tagged_normal_bam = "analysis/bam/{sample_n}.{depth_n}.haplotagged.bam",
        hp_tagged_normal_bai = "analysis/bam/{sample_n}.{depth_n}.haplotagged.bam.bai",
        phased_vcf = "analysis/snvs/clair3/{sample_n}.{depth_n}/phased_merge_output.vcf.gz",
        vntr_bed = config['severus']['vntr']
    output:
        vcf = "analysis/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/somatic_SVs/severus_somatic.vcf",
    params:
        outdir = "analysis/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}"
    threads: 20
    resources:
        mem = 64,
        walltime = 48
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/severus
        set -eu 
        severus --target-bam {input.hp_tagged_tumour_bam} \
            --control-bam {input.hp_tagged_normal_bam} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --vntr-bed {input.vntr_bed} \
            --phasing-vcf {input.phased_vcf}
        """


rule index_haplotagged_bam_a:
    input:
        bam = "analysis/bam/{sample}.{depth}.haplotagged.bam",
    output:
        temp("analysis/bam/{sample}.{depth}.haplotagged.bam.bai")
    threads: 8
    resources:
        mem = 10,
        walltime = 8
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools index -@{threads} {input}
        """

rule haplotagging_whatshap_a:
    input:
        genome = config['reference']['file'],
        vcf = lambda wildcards:  "analysis/snvs/clair3/COLO829_BL.{depth}/phased_merge_output.vcf.gz" if wildcards.sample.startswith("COLO829") else \
        "analysis/snvs/clair3/HCC1937_BL.{depth}/phased_merge_output.vcf.gz",
        bam = "analysis/bam/{sample}.{depth}.bam",
    output:
        bam = temp("analysis/bam/{sample}.{depth}.haplotagged.bam"),
    threads: 8
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {Clair3_sif} \
        whatshap haplotag -o {output.bam} \
        --reference {input.genome} {input.vcf} {input.bam} \
        --ignore-read-groups --tag-supplementary --skip-missing-contigs \
        --output-threads={threads}
        """

rule call_germline_snv_clair3_a:
    input:
        bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        bai = "analysis/bam/{sample_n}.{depth_n}.bam.bai",
        genome = config['reference']['file'],
    output:
        "analysis/snvs/clair3/{sample_n}.{depth_n}/phased_merge_output.vcf.gz"
    params:
        outdir = "analysis/snvs/clair3/{sample_n}.{depth_n}",
        model = "/mnt/backedup/home/jiaZ/working/data/ont_models/clair3_models/r1041_e82_400bps_sup_v430"
    threads: 24
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {Clair3_sif} /opt/bin/run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.genome} \
            --sample_name={wildcards.sample_n} \
            --threads={threads} \
            --platform="ont" \
            --model_path={params.model} \
            --enable_phasing \
            --include_all_ctgs \
            --remove_intermediate_dir \
            --output={params.outdir} 
        """


