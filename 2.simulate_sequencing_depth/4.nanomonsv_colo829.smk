### use nanomonsv to call somatic SVs. 
## set --single_bnd --use_racon, so it can identify single breakend sv as well.
## use --qv20 preset parameters for nanopore Q20 chemistry. 
## set --min_indel_size 10 (default 50) to see if it can get smallers sv.maybe not worth it.
## nanomonsv_filter_simple_repeat_svtype: Post filtering of simple repeat and add svtype to the txt results
## nanomonsv_postprocess_sbnd: annotate svs and classify sbnd.
###########

configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
colo829_samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
hcc1937_samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
wildcard_constraints:
    sample = "|".join(colo829_samples + hcc1937_samples),
    depth = "|".join(["60x","45x","30x","15x"])

include: "helper.smk"

rule all:
    input:
        #expand("analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.txt", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.txt", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),
        #expand("analysis/svs/nanomonsv/{sample_t}.60x.COLO829_BL.60x/{sample_t}.60x.COLO829_BL.60x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #expand("analysis/svs/nanomonsv/{sample_t}.45x.COLO829_BL.45x/{sample_t}.45x.COLO829_BL.45x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #expand("analysis/svs/nanomonsv/{sample_t}.30x.COLO829_BL.30x/{sample_t}.30x.COLO829_BL.30x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #expand("analysis/svs/nanomonsv/{sample_t}.15x.COLO829_BL.15x/{sample_t}.15x.COLO829_BL.15x.nanomonsv.result.txt", sample_t=colo829_samples[:-1])
        #
        #expand("analysis/svs/nanomonsv/{sample_t}.60x.COLO829_BL.45x/{sample_t}.60x.COLO829_BL.45x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #expand("analysis/svs/nanomonsv/{sample_t}.45x.COLO829_BL.30x/{sample_t}.45x.COLO829_BL.30x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #expand("analysis/svs/nanomonsv/{sample_t}.30x.COLO829_BL.15x/{sample_t}.30x.COLO829_BL.15x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        #
        expand("analysis/svs/nanomonsv/{sample_t}.60x.COLO829_BL.30x/{sample_t}.60x.COLO829_BL.30x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),
        expand("analysis/svs/nanomonsv/{sample_t}.45x.COLO829_BL.15x/{sample_t}.45x.COLO829_BL.15x.nanomonsv.result.txt", sample_t=colo829_samples[:-1]),

misc = config['nanomonsv']['misc_scripts_path']

rule call_somatic_sv_nanomonsv_get_a:
    input:
        [f"analysis/svs/nanomonsv/{{sample_t}}.{{depth_t}}/{{sample_t}}.{{depth_t}}.{type}.sorted.bed.gz" for type in ['bp_info','deletion','insertion']],
        "analysis/svs/nanomonsv/{sample_t}.{depth_t}/{sample_t}.{depth_t}.rearrangement.sorted.bedpe.gz",
        [f"analysis/svs/nanomonsv/{{sample_n}}.{{depth_n}}/{{sample_n}}.{{depth_n}}.{type}.sorted.bed.gz" for type in ['bp_info','deletion','insertion']],
        "analysis/svs/nanomonsv/{sample_n}.{depth_n}/{sample_n}.{depth_n}.rearrangement.sorted.bedpe.gz",
        genome = config['reference']['file'],
        tumour_bam = "analysis/bam/{sample_t}.{depth_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.{depth_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.{depth_n}.bam",
        control_panel_path = config['nanomonsv']['control_panel_path'],
    output:
        txt = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.txt",
        vcf = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.vcf",
        sbnd_txt = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.sbnd.result.txt",
        sread = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.supporting_read.txt",
    params:
        tumour_prefix = "analysis/svs/nanomonsv/{sample_t}.{depth_t}/{sample_t}.{depth_t}",
        normal_prefix = "analysis/svs/nanomonsv/{sample_n}.{depth_n}/{sample_n}.{depth_n}",
        panel_prefix="hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control",
        final_outdir = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}"
    threads: 20
    resources:
        mem = 48,
        walltime = 36
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        nanomonsv get {params.tumour_prefix} {input.tumour_bam} {input.genome} \
        --control_prefix {params.normal_prefix} --control_bam {input.normal_bam} \
        --single_bnd --use_racon --min_indel_size 10 --qv20 \
        --control_panel_prefix {input.control_panel_path}/{params.panel_prefix} --processes {threads}

        mkdir -p {params.final_outdir}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}.{wildcards.depth_t}/{wildcards.sample_t}.{wildcards.depth_t}.nanomonsv.result.txt {output.txt}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}.{wildcards.depth_t}/{wildcards.sample_t}.{wildcards.depth_t}.nanomonsv.result.vcf {output.vcf}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}.{wildcards.depth_t}/{wildcards.sample_t}.{wildcards.depth_t}.nanomonsv.sbnd.result.txt {output.sbnd_txt}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}.{wildcards.depth_t}/{wildcards.sample_t}.{wildcards.depth_t}.nanomonsv.supporting_read.txt {output.sread}
        """

rule nanomonsv_parse_a:
    input:
        bam = "analysis/bam/{sample}.{depth}.bam",
        bai="analysis/bam/{sample}.{depth}.bam.bai",
        #genome = config['reference']['file']   # nnmsv supports cram from v0.7.0. for cram input, --reference_fasta is recommended.
    output:
        "analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}.bp_info.sorted.bed.gz",
        "analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}.bp_info.sorted.bed.gz.tbi",
        "analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}.deletion.sorted.bed.gz",
        "analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}.insertion.sorted.bed.gz",
        "analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}.rearrangement.sorted.bedpe.gz",
    params:
        prefix="analysis/svs/nanomonsv/{sample}.{depth}/{sample}.{depth}"
    threads: 1
    resources:
        mem = 24,
        walltime = 48
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        nanomonsv parse {input.bam} {params.prefix}
        """
