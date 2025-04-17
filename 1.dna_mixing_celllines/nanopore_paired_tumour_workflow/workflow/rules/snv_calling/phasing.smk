Clair3_sif = config['clair3']['sif']

rule haplotagging_whatshap:
    input:
        genome = config['reference']['file'],
        ## vcf = "analysis/snvs/clair3/{flowcell}/{mode}/{sample}/phased_merge_output.vcf.gz",
        vcf = lambda wildcards: get_phased_vcf(wildcards.sample, wildcards.flowcell, wildcards.mode),
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        #bam="analysis/snvs/clair3/{flowcell}/{mode}/{sample}.haplotagged.bam"
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.haplotagged.bam"
    log:
        "logs/whatshap/{flowcell}.{mode}.{sample}.haplotagging_whatshap.log"
    benchmark:
        "benchmarks/whatshap/{flowcell}.{mode}.{sample}.haplotagging_whatshap.benchmark.txt"
    threads: 8
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {Clair3_sif} \
        whatshap haplotag -o {output.bam} \
        --reference {input.genome} {input.vcf} {input.bam} \
        --ignore-read-groups --tag-supplementary --skip-missing-contigs \
        --output-threads={threads} &> {log}
        """

rule index_haplotagged_bam:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.haplotagged.bam",
    output:
        "analysis/bam/{flowcell}/{mode}/{sample}.haplotagged.bam.bai",
    threads: 8
    resources:
        mem = 10,
        walltime = 8
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools index -@8 {input}
        """