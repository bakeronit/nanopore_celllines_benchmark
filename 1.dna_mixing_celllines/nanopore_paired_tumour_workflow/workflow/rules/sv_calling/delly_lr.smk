## use Delly to call somatic SV from ONT
## This workflow only take one pair of tumour-normal samples as input 
## you might consider building a multi-sample control panel with all normal samples if you have a cohort.
delly = config['delly']['path']

rule call_somatic_sv_delly:  # paired-samples analysis
    input:
        tumour_bam="analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumour_bai="analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai="analysis/bam/{flowcell}/{mode}/{sample_n}.bam.bai",
        genome=config['reference']['file']
    params:
        excl = f"-x {config['delly']['bed']}" if config['delly']['bed'] != None else ""
    threads: 2
    resources:
        mem = 48,
        walltime = 48
    output:
        bcf = "analysis/svs/delly/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.pre.bcf",
        samples_tsv = "analysis/svs/delly/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.tsv",
        filtered_vcf = "analysis/svs/delly/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.vcf"
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
            {output.bcf} | bcftools view -Ov > {output.filtered_vcf}
        """