#DP_somatic_sif = "/mnt/backedup/home/jiaZ/working/imgs/deepvariant_head561058313.sif"
DP_somatic_sif = "/mnt/backedup/home/jiaZ/working/imgs/deepsomatic/deepsomatic_1.6.0.sif"

rule call_somatic_snv_deepsomatic:
    input:
        genome=config['reference']['file'],
        tumor_bam = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumor_bai = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai"
    output:
        "analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.vcf.gz"
    benchmark:
        "benchmarks/deepsomatic/{flowcell}.{mode}.{sample_t}.{sample_n}.benchmark.txt"
    params:
        model="/mnt/backedup/home/jiaZ/working/data/ont_models/dpsomatic_model/weights-143-0.987994.ckpt"
    threads: 48
    envmodules:
        "singularity/3.7.1"
    resources:
        mem=40,
        walltime=100
    shell:
        """
        singularity exec {DP_somatic_sif} /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --ref={input.genome} \
        --model_type=ONT_R104 \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --num_shards={threads} \
        --logging_dir=logs/deepsomatic/{wildcards.flowcell}.{wildcards.mode}.{wildcards.sample_t}.{wildcards.sample_n} \
        --customized_model={params.model}
        """
