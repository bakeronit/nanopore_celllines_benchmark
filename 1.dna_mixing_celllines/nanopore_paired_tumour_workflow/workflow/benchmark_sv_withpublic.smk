workdir: "../../work"

configfile: "../config/config.yaml"

include: "rules/common.smk"
pairs = generate_paired_samples(samples_df)

rule all:
    input:
        [f"analysis/benchmark_public/svs/nanomonsv/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark_public/svs/severus/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark_public/svs/savana/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark_public/svs/delly/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],


rule jasmine_merge_a:
    input:
        ref = config['reference']['file'],
        vcf = "analysis/svs/simple/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.vcf",
    output:
        "analysis/benchmark_public/svs/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.merged.vcf",
    params:
        gs = "colo829.refined.vcf",
        tmp_dir = lambda wildcards: f"analysis/benchmark_public/svs/tmp/{wildcards.tool}",
        file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.file.txt",
        refined_vcf = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.refined.vcf",
        refined_file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.refined_file.txt",
        out_dir = lambda wildcards: f"analysis/benchmark_public/svs/{wildcards.tool}/{wildcards.flowcell}/{wildcards.mode}"
    threads: 1
    resources:
        mem = 1,
        walltime = 1
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/jasmine
        set -eu

        mkdir -p {params.tmp_dir}
        ls {input.vcf} > {params.tmp_dir}/{params.file_list}

        jasmine --preprocess_only --pre_normalize --dup_to_ins file_list={params.tmp_dir}/{params.file_list} genome_file={input.ref} out_dir={params.tmp_dir}
        jasmine file_list={params.tmp_dir}/{wildcards.sample_t}.{wildcards.sample_n}_dupToIns_normalizeTypes.vcf out_dir={params.tmp_dir} \
        --comma_filelist max_dist=200 --use_end --ignore_strand --allow_intrasample --nonlinear_dist out_file={params.tmp_dir}/{params.refined_vcf}

        ls -U {params.gs} {params.tmp_dir}/{params.refined_vcf} > {params.tmp_dir}/{params.refined_file_list}
        jasmine file_list={params.tmp_dir}/{params.refined_file_list} out_file={output} --ignore_strand --use_end

        jasmine --dup_to_ins --postprocess_only out_file={output} out_dir={params.out_dir}
        """
