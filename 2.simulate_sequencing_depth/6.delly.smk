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
delly = config['delly']['path']

rule all:
    input:
        expand("analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf", my_dirty_combinator,sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf", my_dirty_combinator,sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"]),

rule call_somatic_sv_delly_a:
    input:
        genome = config['reference']['file'],
        tumour_bam = "analysis/bam/{sample_t}.{depth_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.{depth_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.{depth_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.{depth_n}.bam.bai",
    output:
        bcf = "analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.pre.bcf",
        samples_tsv = "analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.tsv",
        vcf = "analysis/svs/delly/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.vcf",
    params:
        excl = f"-x {config['delly']['bed']}" if config['delly']['bed'] != None else ""
    threads: 2
    resources:
        mem = 48,
        walltime = 48
    envmodules:
        "samtools/1.17",
        "bcftools/1.16"
    shell:
        """
        export OMP_NUM_THREADS={threads} # delly primarily parallises on the sample level.
        {delly} lr \
            -t ALL \
            -y ont \
            -o {output.bcf} {params.excl} \
            -g {input.genome} {input.tumour_bam} {input.normal_bam}
        
        function get_sample_id() {{
            echo "$(samtools view -H ${{1}} | perl -lne 'print ${{1}} if /\\sSM:(\\S+)/' | head -n 1 )"
        }}

        TID=$(get_sample_id "{input.tumour_bam}")
        CID=$(get_sample_id "{input.normal_bam}")
        printf "${{TID}}\\ttumor\\n${{CID}}\\tcontrol\\n" > {output.samples_tsv}

        {delly} filter \
            -f somatic \
            -p \
            -s {output.samples_tsv} \
            {output.bcf} | bcftools view -Ov > {output.vcf}
        """
