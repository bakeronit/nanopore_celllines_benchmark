##### setup report #####


#samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
include: "helper.smk"


rule all:
    input:
        #expand("analysis/benchmark/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['COLO829_BL'], depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/benchmark/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['COLO829_BL'], depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/benchmark/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['COLO829_BL'], depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/benchmark/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['COLO829_BL'], depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['HCC1937_BL'], depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['HCC1937_BL'], depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['HCC1937_BL'], depth_n=["60x","45x","30x","15x"]),
        expand("analysis/benchmark/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf", my_dirty_combinator, sample_t=samples[:-1], depth_t=["60x","45x","30x","15x"], sample_n=['HCC1937_BL'], depth_n=["60x","45x","30x","15x"]),


def get_somatic_sv(wildcards):
    match wildcards.tool:
        case 'nanomonsv':
            return "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.pass.svtype.txt"
        case 'severus':
            return "analysis/svs/severus/{sample_t}.{depth_t}.{sample_n}.{depth_n}/somatic_SVs/severus_somatic.vcf"
        case 'delly':
            return "analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf"
        case 'savana':
            return "analysis/svs/savana/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.classified.somatic.vcf"
        case _:
            raise VauleError("tool can only be nanomonsv, severus, savana, and delly") 

rule get_simple_svtype:
    input:
        get_somatic_sv
    output:
       simple = "analysis/svs/simple/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.rawlen50.vcf",
       filtered = "analysis/svs/simple/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf"
    threads: 1
    params:
        script = "~/working/dev/SV_vcf/simple_event_annotation.py",
        ref = lambda w: "-r ~/working/data/genome/reference.fasta" if w.tool=="nanomonsv" else ""
    resources:
        mem = 1,
        walltime = 1
    envmodules:
        "conda-envs/base",
        "bcftools/1.19"
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/jasmine
        set -eu

        python {params.script} {input} -t {wildcards.tool} -o {output.simple} {params.ref}

        bcftools filter -i '((SVTYPE="DUP" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="INS") && SVLEN>50) || (SVTYPE="TRA")' {output.simple} | grep -Ew "^#|^chr([1-9]|1[0-9]|2[0-2]|X|Y)" > {output.filtered}
        """

rule jasmine_merge:
    input:
        ref = "/mnt/backedup/home/jiaZ/working/data/genome/reference.fasta",
        vcf = "analysis/svs/simple/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf"
    threads: 1
    output:
        "analysis/benchmark/svs/{tool}/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.merged.vcf"
    params:
        gs = lambda wildcards: f"/mnt/backedup/home/jiaZ/working/general/goldstandard/structural_variation/analysis/svs/jasmine_merge/{wildcards.sample_t.lower().split("_")[0]}/{wildcards.sample_t.lower().split("_")[0]}.final_merged.supp2_dupToIns.vcf",
        tmp_dir = lambda wildcards: f"analysis/benchmark/svs/tmp/{wildcards.tool}",
        file_list = lambda wildcards: f"{wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n}.file.txt",
        refined_vcf = lambda wildcards: f"{wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n}.refined.vcf",
        refined_file_list = lambda wildcards: f"{wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n}.refined_file.txt",
        out_dir = lambda wildcards: f"analysis/benchmark/svs/{wildcards.tool}/{wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n}"
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
        jasmine file_list={params.tmp_dir}/{wildcards.sample_t}.{wildcards.depth_t}.{wildcards.sample_n}.{wildcards.depth_n}_dupToIns_normalizeTypes.vcf out_dir={params.tmp_dir} \
        --comma_filelist max_dist=200 --use_end --ignore_strand --allow_intrasample --nonlinear_dist out_file={params.tmp_dir}/{params.refined_vcf}

        ls -U {params.gs} {params.tmp_dir}/{params.refined_vcf} > {params.tmp_dir}/{params.refined_file_list}
        /mnt/backedup/home/jiaZ/working/local/Jasmine/jasmine file_list={params.tmp_dir}/{params.refined_file_list} out_file={output} out_dir={params.out_dir} --ignore_strand --use_end 

        jasmine --dup_to_ins --postprocess_only out_file={output} out_dir={params.out_dir}
        """
