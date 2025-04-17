script="/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/scripts"
bamcov = config['bamcov']
rule bc_summary:
    input:
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam"
    output:
        summary = "analysis/qc/basecalling/{flowcell}/{mode}/{sample}/{run}.summary.txt",
        dy = "analysis/qc/basecalling/{flowcell}/{mode}/{sample}/{run}.data_yield.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 4
    shell:
        """
        {DORADO} summary {input} > {output.summary}
        {script}/data_yield.awk {output.summary} > {output.dy}
        """

rule bam_stats:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.mosdepth.global.dist.txt",
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.mosdepth.summary.txt",
    params:
        prefix = "analysis/qc/bam/{flowcell}/{mode}/{sample}"
    threads: 4
    resources:
        mem = 10,
        walltime = 8
    envmodules:
        "mosdepth/0.2.9"
    shell:
        """
        mosdepth -n -t {threads} {params.prefix} {input}
        """

rule bam_cov:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.bamcov.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 12
    shell:
        """
        {bamcov} -H {input} -o {output}
        """

rule align_N50:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.bamN50.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 2
    envmodules:
        "python/3.9.13"
    shell:
        """
        python {script}/getAlignmentN50.py {input} > {output}
        """