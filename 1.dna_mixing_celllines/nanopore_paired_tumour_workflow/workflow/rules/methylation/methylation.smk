MODKIT = config['modkit']

rule bamTobedmethyl:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai",
        genome = config['reference']['file'],
    output:
        gz="analysis/mod/{flowcell}/{mode}/{sample}.bed.gz",
        tbi="analysis/mod/{flowcell}/{mode}/{sample}.bed.gz.tbi"
    params:
        bed = "analysis/mod/{flowcell}/{mode}/{sample}.bed"
    threads: 12
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{flowcell}.{mode}.{sample}.pileup.log"
    benchmark:
        "benchmarks/modkit/{flowcell}.{mode}.{sample}.benchmark.txt"
    envmodules:
        "htslib/1.17"
    shell:
        """
        {MODKIT} pileup {input.bam} {params.bed} --ref {input.genome} \
        --threads {threads} --combine-strands --cpg --log-filepath {log}
        bgzip {params.bed}

        tabix -p bed {output.gz}
        """

rule modsummary:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    output:
        "analysis/mod/{flowcell}/{mode}/{sample}.mod_summary.txt"
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{flowcell}.{mode}.{sample}.summary.log"
    shell:
        """
        {MODKIT} summary -t {threads} {input.bam} --log-filepath {log} > {output}
        """