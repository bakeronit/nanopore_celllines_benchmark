Clair3_sif = config['clair3']['sif']
DEEPVARIANT_sif = config['deepvariant']
ClairS_sif = config['clairS']['sif']
PEPPER_sif = config['pepper']['sif']

DP_somatic_sif = "/mnt/backedup/home/jiaZ/working/imgs/deepvariant_head561058313.sif"

rule clair3_snv_calling:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        csi = "analysis/bam/{flowcell}/{mode}/{sample}.bam.csi",
        genome = config['reference']['file'],
    output:
        "analysis/snvs/clair3/{flowcell}/{mode}/{sample}/merge_output.vcf.gz"
    params:
        outdir = "analysis/snvs/clair3/{flowcell}/{mode}/{sample}",
        model=get_clair3_model
    log:
        "logs/clair3/{flowcell}.{mode}.{sample}.log"
    benchmark:
        "benchmarks/clair3/{flowcell}.{mode}.{sample}.benchmark.txt"
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
            --sample_name={wildcards.sample} \
            --ref_fn={input.genome} \
            --threads={threads} \
            --platform="ont" \
            --model_path={params.model} \
            --enable_phasing \
            --output={params.outdir} 2>{log}
        """

## deepvariant only available for R10
rule call_germline_snp_deepvariant:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        vcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.vcf.gz",
        gvcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.g.vcf.gz"
    log:
        "logs/deepvariant/{flowcell}.{mode}.{sample}.log"
    benchmark:
        "benchmarks/deepvariant/{flowcell}.{mode}.{sample}.benchmark.txt"
    threads: 12
    resources:
        gpu = 1,
        mem = 48,
        walltime = 48
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        singularity run --nv {DEEPVARIANT_sif} \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type ONT_R104 \
            --ref {input.genome} \
            --reads {input.bam} \
            --sample_name {wildcards.sample} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {threads}  2>{log}
        """

rule call_germline_snp_pepper:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        vcf="analysis/snvs/pepper/{flowcell}/{mode}/{sample}/{sample}.vcf.gz"
    log:
        "logs/pepper/{flowcell}.{mode}.{sample}.log"
    benchmark:
        "benchmarks/pepper/{flowcell}.{mode}.{sample}.benchmark.txt"
    threads: 24
    resources:
        mem = 48,
        walltime = 48
    envmodules:
        "singularity/3.7.1"
    params:
        output_dir = "analysis/snvs/pepper/{flowcell}/{mode}/{sample}",
        mode = lambda w: '--ont_r9_guppy5_sup' if w.flowcell == 'R9' else '--ont_r10_q20'
    shell:
        """
        singularity exec {PEPPER_sif} \
        run_pepper_margin_deepvariant call_variant \
            --fasta {input.genome} \
            --bam {input.bam} \
            --output_dir {params.output_dir} \
            --output_prefix {wildcards.sample} \
            --sample_name {wildcards.sample} \
            --phased_output --skip_final_phased_bam \
            {params.mode} \
            --threads {threads}  2>{log}
        """


### run clairS with paired tumor-normal bam files
rule clairS_snv_calling:
    input:
        genome = config['reference']['file'],
        tumor_bam = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumor_csi = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam.csi",
        normal_bam = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_csi = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam.csi"
    output:
        "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/output.vcf.gz"
    params:
        platform = lambda w: config['clairS']['platform'][w.flowcell],
        outdir = "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}",
        indel_option = lambda w: "--enable_indel_calling" if w.flowcell == 'R10' else ""
    log:
        "logs/clairS/{flowcell}.{mode}.{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/clairS/{flowcell}.{mode}.{sample_t}.{sample_n}.benchmark.txt"
    threads: 24
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {ClairS_sif} /opt/bin/run_clairs \
            --tumor_bam_fn {input.tumor_bam} \
            --normal_bam_fn {input.normal_bam} \
            --ref_fn {input.genome} \
            --include_all_ctgs {params.indel_option}\
            --sample_name {wildcards.sample_t} \
            --threads {threads} \
            --platform {params.platform} \
            --output_dir {params.outdir} \
            --remove_intermediate_dir \
            --conda_prefix /opt/conda/envs/clairs 2>{log}
        """

rule deepSomatic:
    input:
        genome=config['reference']['file'],
        tumor_bam = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumor_csi = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam.csi",
        normal_bam = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_csi = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam.csi"
    output:
        "analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.vcf.gz"
    log:
        "logs/deepsomatic/{flowcell}.{mode}.{sample_t}.{sample_n}.log"
    params:
        model="/mnt/backedup/home/jiaZ/working/data/ont_models/dpsomatic_model/weights-143-0.987994.ckpt"
    threads: 12
    envmodules:
        "singularity/3.7.1"
    resources:
        mem=30,
        walltime=48
    shell:
        """
        singularity exec {DP_somatic_sif} /opt/deepvariant/bin/run_deepsomatic \
        --model_type=ONT_R104 \
        --ref={input.genome} \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --use_keras_model \
        --num_shards={threads} \
        --logging_dir=logs/deepsomatic \
        --customized_model={params.model}
        """   
