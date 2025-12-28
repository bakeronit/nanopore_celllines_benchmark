#!/usr/bin/env python3

rule delly_call:
    input:
        ref = config['reference'],
        tumor_bam = "analysis/bam/{cell}/{lib}/{tumor}.bam",
        normal_bam = "analysis/bam/{cell}/{lib}/{normal}.bam",
        blacklist_bed = config['blacklist_bed'],
    output:
        bcf = "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.pre.bcf",
        sample_tsv = "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.tsv",
    params:
        delly_path = config['delly_path']
    threads: 2
    resources:
        mem = 30,
        walltime = 48
    envmodules:
        "bcftools/1.19",
        "samtools/1.17"
    shell:
        """
        export OMP_NUM_THREADS={threads}

        {params.delly_path}/delly call -x {input.blacklist_bed} -t ALL \
        -o {output.bcf} -g {input.ref} \
            {input.tumor_bam} {input.normal_bam}

        function get_sample_id() {{
            echo "$(samtools view -H ${{1}} | perl -lne 'print ${{1}} if /\\sSM:(\\S+)/' | head -n 1 )"
        }}

        TID=$(get_sample_id "{input.tumor_bam}")
        CID=$(get_sample_id "{input.normal_bam}")

        printf "${{TID}}\\ttumor\\n${{CID}}\\tcontrol\\n" > {output.sample_tsv}
        """

rule delly_filter:
    input:
        bcf = "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.pre.bcf",
        sample_tsv = "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.tsv"
    output:
        "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.vcf"
    params:
        delly_path = config['delly_path']
    threads: 1
    resources:
        mem = 5,
        walltime = 10
    envmodules:
        "bcftools/1.19"
    shell:
        """
        {params.delly_path}/delly filter -f somatic -s {input.sample_tsv} {input.bcf} | bcftools view -Ov > {output}
        """

