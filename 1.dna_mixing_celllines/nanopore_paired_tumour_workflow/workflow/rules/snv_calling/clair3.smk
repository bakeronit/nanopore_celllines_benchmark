## germline calling with clair3
## --enable_phasing will generates phased vcf and haplotagged bam files using whatshap.
## --include_all_ctgs to call variants for non-standard chr as well.
## --remove_intermediate_dir will remove intermediate directory including bam, etc. /tmp
## --enable_long_indel will call long indels (>50bp) consider using it in the future.
## after calling, we used whatshap to haplotag bam file with phased vcf
#######

Clair3_sif = config['clair3']['sif']

rule call_germline_snv_clair3:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai",
        genome = config['reference']['file'],
    output:
        "analysis/snvs/clair3/{flowcell}/{mode}/{sample}/phased_merge_output.vcf.gz"
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
            --ref_fn={input.genome} \
            --sample_name={wildcards.sample} \
            --threads={threads} \
            --platform="ont" \
            --model_path={params.model} \
            --enable_phasing \
            --include_all_ctgs \
            --remove_intermediate_dir \
            --output={params.outdir} | tee -a {log}
        """