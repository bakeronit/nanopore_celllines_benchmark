#!/usr/bin/env python3

rule gridss:
    input:
        ref = config['reference'],
        tumor_bam = "analysis/bam/{cell}/{lib}/{tumor}.bam",
        normal_bam = "analysis/bam/{cell}/{lib}/{normal}.bam",
        blacklist_bed = "/mnt/backedup/home/jiaZ/working/data/hg38-blacklist.v2.bed"
    output:
        raw_vcf = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}.raw.vcf",
        filtered_vcf = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}.pass.vcf"
    params:
        workdir = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}",
        gridss_path = config['gridss']['path'],
        gridss_dependencies = config['gridss']['dependencies']
    threads: 24
    resources:
        mem = 30,
        walltime = 48
    envmodules:
        "bcftools/1.19",
        "samtools/1.17",
        "R/4.3.1",
        "bwa/0.7.15",
        "kraken2/2.1.2",
        "RepeatMasker/4.1.0"
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        ln -s {input.ref} {params.workdir}/reference.fasta
        
        bwa index {params.workdir}/reference.fasta

        {params.gridss_path}/gridss -j {params.gridss_dependencies} \
            --reference {params.workdir}/reference.fasta \
            --threads {threads} \
            --workingdir {params.workdir} \
            --blacklist {input.blacklist_bed} \
            --assembly {params.workdir}/gridss_assembly.bam \
            --output {output.raw_vcf} \
            {input.normal_bam} \
            {input.tumor_bam}

        ## get somatic SV + filtering
        bcftools filter -i "FORMAT/QUAL[0] == 0 && FILTER == '.'" \
            -Ov -o {output.filtered_vcf} {output.raw_vcf} 
        """

rule gripss_filter:
    input:
        ref = config['reference'],
        vcf = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}.raw.vcf",
    output:
        "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}/{tumor}.bam.gripss.filtered.vcf.gz"
    params:
        gridss_path = config['gridss']['path'],
        pon_sgl = "/mnt/backedup/home/jiaZ/working/local/gridss/v5_34/ref/38/sv/sgl_pon.38.bed.gz",
        pon_sv = "/mnt/backedup/home/jiaZ/working/local/gridss/v5_34/ref/38/sv/sv_pon.38.bedpe.gz",
        repeat_mask = "/mnt/backedup/home/jiaZ/working/local/gridss/v5_34/ref/38/sv/repeat_mask_data.38.fa.gz",
        outdir = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}"
    threads: 1
    resources:
        mem = 10,
        walltime = 1
    shell:
        """
        java -jar {params.gridss_path}/gripss_v2.4.jar \
            -sample {wildcards.tumor}.bam -reference {wildcards.normal}.bam -ref_genome_version 38 \
            -ref_genome {input.ref} \
            -vcf {input.vcf} \
            -pon_sgl_file {params.pon_sgl} \
            -pon_sv_file {params.pon_sv}  \
            -repeat_mask_file {params.repeat_mask} \
            -output_dir {params.outdir}
        """
