## germline variant calling with pepper-margin-deepvariant pipeline
## it will generates phased vcf and haplotagged bam files
## set --phased_output to enable phasing
## set --skip_final_phased_bam will remove the bam files
##### 

PEPPER_MARGIN_DV_sif = config['pepper']['sif']

rule call_germline_snv_pepper:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai="analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    output:
        vcf="analysis/snvs/pepper/{flowcell}/{mode}/{sample}/{sample}.phased.vcf.gz"
    log:
        "logs/pepper/{flowcell}.{mode}.{sample}.log"
    benchmark:
        "benchmarks/pepper/{flowcell}.{mode}.{sample}.benchmark.txt"
    threads: 24
    resources:
        mem = 60,
        walltime = 48
    envmodules:
        "singularity/3.7.1"
    params:
        output_dir = "analysis/snvs/pepper/{flowcell}/{mode}/{sample}",
        mode = lambda w: '--ont_r9_guppy5_sup' if w.flowcell == 'R9' else '--ont_r10_q20'
    shell:
        """
        singularity exec {PEPPER_MARGIN_DV_sif} \
        run_pepper_margin_deepvariant call_variant \
            --bam {input.bam} \
            --fasta {input.genome} \
            --output_dir {params.output_dir} \
            --output_prefix {wildcards.sample} \
            --sample_name {wildcards.sample} \
            --phased_output \
            {params.mode} \
            --threads {threads} | tee -a {log}
        """