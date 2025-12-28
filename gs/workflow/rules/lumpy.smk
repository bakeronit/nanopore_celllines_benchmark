#!/usr/bin/env python3

rule get_discordants:
    input:
        bam = "analysis/bam/{cell}/{lib}/{sample}.bam"    
    output:
        bam = "analysis/bam/{cell}/{lib}/{sample}.discordants.bam",
    threads: 10
    resources:
        mem = 20,
        walltime = 10
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools view -@{threads} -b -F 1294 {input.bam} |samtools sort -@{threads} > {output.bam}
        """

rule get_splitters:
    input:
        bam = "analysis/bam/{cell}/{lib}/{sample}.bam"
    output:
        bam = "analysis/bam/{cell}/{lib}/{sample}.splitters.bam"
    params:
        lumpy_scripts = config['lumpy']['scripts']
    threads: 10
    resources:
        mem = 20,
        walltime = 10
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools view -@{threads} -h {input.bam} | {params.lumpy_scripts}/extractSplitReads_BwaMem -i stdin |\
        samtools view -Sb - |samtools sort -@{threads} > {output.bam}
        """

rule lumpy:
    input:
        ref = config['reference'],
        tumor_bam = "analysis/bam/{cell}/{lib}/{tumor}.bam",
        normal_bam = "analysis/bam/{cell}/{lib}/{normal}.bam",
        tumor_ds_bam = "analysis/bam/{cell}/{lib}/{tumor}.discordants.bam",
        normal_ds_bam = "analysis/bam/{cell}/{lib}/{normal}.discordants.bam",
        tumor_ss_bam = "analysis/bam/{cell}/{lib}/{tumor}.splitters.bam",
        normal_ss_bam = "analysis/bam/{cell}/{lib}/{normal}.splitters.bam",
        blacklist_bed = config['backlist_bed'],
    output:
        vcf = "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.raw.vcf",
    params:
        workdir = "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}",
        lumpy_path = config['lumpy']['bin'],
    threads: 1
    resources:
        mem = 30,
        walltime = 48
    envmodules:
        "python/2.7.13",
        "bcftools/1.19",
        "samtools/1.17",
        "sambamba/1.0.1",
        "samblaster/0.1.26",
    shell:
        """
        {params.lumpy_path}/lumpyexpress \
            -B {input.tumor_bam},{input.normal_bam} \
            -D {input.tumor_ds_bam},{input.normal_ds_bam} \
            -S {input.tumor_ss_bam},{input.normal_ss_bam} \
            -x {input.blacklist_bed} \
            -o {output.vcf}
        """

#rule lumpy_filter:
#    input:
#        "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.raw.vcf",
#    output:
#        "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.pass.vcf",
#    threads: 1
#    resources:
#        mem = 2,
#        walltime = 1
#    envmodules:
#        "bcftools/1.19"
#    shell:
#        """
#        ## get somatic SV + filtering
#        bcftools filter -i "FORMAT/SU[1] == 0 && FILTER == '.'" \
#            -Ov -o {output} {input}
#        """

rule lumpy_svtyper:
    input:
        vcf = "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.raw.vcf",
        tumor_bam = "analysis/bam/{cell}/{lib}/{tumor}.bam",
        normal_bam = "analysis/bam/{cell}/{lib}/{normal}.bam",
    output:
        vcf = "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.gt.vcf",
    threads: 1
    resources:
        mem = 10,
        walltime = 24
    envmodules:
        "python/2.7.13"
    shell:
        """
        svtyper \
         -i {input.vcf} -B {input.tumor_bam},{input.normal_bam}  > {output.vcf} 
        """

rule lumpy_filter:
    input:
        "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.gt.vcf"
    output:
        "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.somatic.vcf"
    threads: 1
    resources:
        mem = 1,
        walltime = 1
    envmodules:
        "bcftools/1.19"
    shell:
        """
        bcftools filter -i 'GT[0] ~ "1"  && AO[1:*] = 0' {input} > {output}
        """
