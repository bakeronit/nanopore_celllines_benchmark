### run clairS with paired tumor-normal bam files
ClairS_sif = config['clairS']['sif']
rule call_somatic_snv_clairS:
    input:
        genome = config['reference']['file'],
        tumor_bam = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumor_bai = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai"
    output:
        "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/output.vcf.gz"
    params:
        platform = lambda w: config['clairS']['platform'][w.flowcell],
        outdir = "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}",
        clair3_model = lambda w: config['clair3']['model_r10'] if w.flowcell == 'R10' else config['clair3']['model_r9'],
        indel_option = lambda w: "--enable_indel_calling" if w.flowcell == 'R10' else ""   # indel calling currently only available for R10 data.
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
            --clair3_model_path {params.clair3_model} \
            --include_all_ctgs {params.indel_option} \
            --sample_name {wildcards.sample_t} \
            --threads {threads} \
            --platform {params.platform} \
            --output_dir {params.outdir} \
            --remove_intermediate_dir \
            --conda_prefix /opt/conda/envs/clairs | tee -a {log}
        """
