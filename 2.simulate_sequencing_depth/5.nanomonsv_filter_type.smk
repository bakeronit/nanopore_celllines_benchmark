configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
colo829_samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
hcc1937_samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]

wildcard_constraints:
    sample = "|".join(colo829_samples + hcc1937_samples),
        depth = "|".join(["60x","45x","30x","15x"])

include: "helper.smk"

rule all:
    input:
        expand("analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.pass.svtype.txt", my_dirty_combinator, sample_t=colo829_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['COLO829_BL'],depth_n=["60x","45x","30x","15x"]),
        expand("analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.pass.svtype.txt", my_dirty_combinator, sample_t=hcc1937_samples[:-1], depth_t=["60x","45x","30x","15x"],sample_n=['HCC1937_BL'],depth_n=["60x","45x","30x","15x"])

misc = config['nanomonsv']['misc_scripts_path']

rule nanomonsv_filter_simple_repeat_svtype:
    input:
        result = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.txt",
        simple_repeat = config['nanomonsv']['simple_repeat']
    output:
        filt = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.txt",
        passed = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.pass.txt",
        svtype = "analysis/svs/nanomonsv/{sample_t}.{depth_t}.{sample_n}.{depth_n}/{sample_t}.{depth_t}.{sample_n}.{depth_n}.nanomonsv.result.filt.pass.svtype.txt"
    threads: 1
    resources:
        mem = 2,
        walltime = 1
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        python3 {misc}/add_simple_repeat.py {input.result} {output.filt} {input.simple_repeat}

        head -n 1 {output.filt} > {output.passed} 
        tail -n +2 {output.filt} |grep PASS >> {output.passed}

        python3 {misc}/sv_type.py {output.passed} {output.svtype}
        """
