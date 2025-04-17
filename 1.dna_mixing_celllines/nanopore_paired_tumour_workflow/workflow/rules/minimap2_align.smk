rule align_minimap2:
    input:
        ubam = "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.trimmed_repaired.ubam" if config['filter']['trim'] else \
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam",
        genome = config['reference']['file']
    output:
        bam = temp("analysis/bam/{flowcell}/{mode}/{sample}/{run}.bam"),
        bai = temp("analysis/bam/{flowcell}/{mode}/{sample}/{run}.bam.bai")
    envmodules:
        "samtools/1.17",
        "minimap2/2.26"
    threads: 24
    resources:
        mem = 36,
        walltime = 24
    params:
        qs = config['filter']['read_qs'],
        rg = "\"@RG\\tID:{sample}.{run}\\tPL:ONT\\tSM:{sample}\""
    benchmark:
        "benchmarks/minimap2/{flowcell}.{sample}.{run}.{mode}.benchmark.txt"
    log:
        "logs/minimap2/{flowcell}.{sample}.{run}.{mode}.log"
    shell: 
        """
        samtools view -e '[qs] >= {params.qs}' {input.ubam} | \
        samtools fastq -@8 -T "*" | \
            minimap2 -R {params.rg} \
            -y --MD -ax map-ont -t {threads} {input.genome} - | \
            samtools sort -@8 -O BAM --write-index -o {output.bam}##idx##{output.bai} - &>{log}
        """ 

rule merge_bam:
    input:
        get_bam_of_runs
    output:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    threads: 10
    resources:
        mem = 20,
        walltime = 8
    run:
        if len(input)>1:
            shell("module load samtools/1.17 && samtools merge -@{threads} --write-index {output.bam}##idx##{output.bai} {input}")
        else:
            shell("ln -sr {input} {output.bam} && ln -sr {input}.bai {output.bai} && touch -h {output.bam} && touch -h {output.bai}")