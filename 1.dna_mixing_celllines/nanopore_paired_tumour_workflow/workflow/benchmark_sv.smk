workdir: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work"

##### setup report #####
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"

include: "rules/common.smk"
include: "rules/sv_calling/nanomonsv.smk"
include: "rules/sv_calling/severus.smk"
include: "rules/sv_calling/savana.smk"
include: "rules/sv_calling/delly_lr.smk"

pairs = generate_paired_samples(samples_df)
wildcard_constraints:
    flowcell = "|".join(['R9','R10']),
    mode = "|".join(['sup','hac'])

rule all:
    input:
        [f"analysis/svs/simple/nanomonsv/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.vcf" for pair in pairs for m in MODE],
        [f"analysis/svs/simple/severus/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.vcf" for pair in pairs for m in MODE],
        [f"analysis/svs/simple/savana/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.vcf" for pair in pairs for m in MODE],
        [f"analysis/svs/simple/delly/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark/svs/nanomonsv/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark/svs/severus/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark/svs/savana/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark/svs/delly/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}.merged.vcf" for pair in pairs for m in MODE],
        [f"analysis/benchmark/svs/{pair['flowcell_version']}_{m}_{pair['tumour']}.{pair['normal']}.all_merged.vcf" for pair in pairs for m in MODE],


def get_somatic_sv(wildcards):
    match wildcards.tool:
        case 'nanomonsv':
            return rules.nanomonsv_filter_simple_repeat_svtype.output.svtype
        case 'severus':
            return rules.call_somatic_sv_severus.output
        case 'delly':
            return rules.call_somatic_sv_delly.output.filtered_vcf
        case 'savana':
            return rules.call_somatic_sv_savana.output[-1]
        case _:
            raise VauleError("tool can only be nanomonsv, severus, savana, and delly") 

rule get_simple_svtype:
    input:
        get_somatic_sv
    output:
       simple = "analysis/svs/simple/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.rawlen50.vcf",
       filtered = "analysis/svs/simple/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.vcf"
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

        bcftools filter -i '((SVTYPE="DUP" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="INS") && SVLEN>=50) || (SVTYPE="TRA")' {output.simple} | grep -Ew "^#|^#CHROM|^chr([1-9]|1[0-9]|2[0-2]|X|Y)" > {output.filtered}
        """

rule jasmine_merge:
    input:
        ref = "/mnt/backedup/home/jiaZ/working/data/genome/reference.fasta",
        vcf = "analysis/svs/simple/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.vcf",
    output:
        "analysis/benchmark/svs/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}.merged.vcf",
    params:
        gs = lambda wildcards: f"/mnt/backedup/home/jiaZ/working/general/goldstandard/structural_variation/analysis/svs/jasmine_merge/{wildcards.sample_t.lower().split("_")[0]}/{wildcards.sample_t.lower().split("_")[0]}.final_merged.supp2_dupToIns.vcf",
        tmp_dir = lambda wildcards: f"analysis/benchmark/svs/tmp/{wildcards.tool}",
        file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.file.txt",
        refined_vcf = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.refined.vcf",
        refined_file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.refined_file.txt",
        out_dir = lambda wildcards: f"analysis/benchmark/svs/{wildcards.tool}/{wildcards.flowcell}/{wildcards.mode}"
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
        /mnt/backedup/home/jiaZ/working/local/Jasmine/jasmine file_list={params.tmp_dir}/{params.refined_file_list} out_file={output} --ignore_strand --use_end

        jasmine --dup_to_ins --postprocess_only out_file={output} out_dir={params.out_dir}
        """

rule jasmine_merge_all:
    input:
        ref = "/mnt/backedup/home/jiaZ/working/data/genome/reference.fasta",
        delly = "analysis/benchmark/svs/tmp/delly/{flowcell}_{mode}_{sample_t}.{sample_n}.refined.vcf",
        nanomonsv = "analysis/benchmark/svs/tmp/nanomonsv/{flowcell}_{mode}_{sample_t}.{sample_n}.refined.vcf",
        savana = "analysis/benchmark/svs/tmp/savana/{flowcell}_{mode}_{sample_t}.{sample_n}.refined.vcf",
        severus = "analysis/benchmark/svs/tmp/severus/{flowcell}_{mode}_{sample_t}.{sample_n}.refined.vcf"
    output:
        tools_merged = "analysis/benchmark/svs/{flowcell}_{mode}_{sample_t}.{sample_n}.tools_merged.vcf",
        supp2 = "analysis/benchmark/svs/{flowcell}_{mode}_{sample_t}.{sample_n}.tools_merged.supp2.vcf",
        all_merged = "analysis/benchmark/svs/{flowcell}_{mode}_{sample_t}.{sample_n}.all_merged.vcf",
    params:
        gs = lambda wildcards: f"/mnt/backedup/home/jiaZ/working/general/goldstandard/structural_variation/analysis/svs/jasmine_merge/{wildcards.sample_t.lower().split("_")[0]}/{wildcards.sample_t.lower().split("_")[0]}.final_merged.supp2_dupToIns.vcf",
        tmp_dir = lambda wildcards: f"analysis/benchmark/svs/tmp",
        file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.all_file.txt",
        bench_file_list = lambda wildcards: f"{wildcards.flowcell}_{wildcards.mode}_{wildcards.sample_t}.{wildcards.sample_n}.all_file.bc.txt",
        out_dir = lambda wildcards: f"analysis/benchmark/svs"
    threads: 1
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

        ls -U {input.delly} {input.nanomonsv} {input.savana} {input.severus} > {params.tmp_dir}/{params.file_list}
        jasmine file_list={params.tmp_dir}/{params.file_list} out_file={output.tools_merged} --ignore_strand --use_end

        bcftools filter -i "SUPP>1" {output.tools_merged} > {output.supp2}

        ls -U {params.gs} {output.supp2} > {params.tmp_dir}/{params.bench_file_list}
        jasmine file_list={params.tmp_dir}/{params.bench_file_list} out_file={output.all_merged} --ignore_strand --use_end
        jasmine --dup_to_ins --postprocess_only out_file={output.all_merged} out_dir={params.out_dir}
        """
