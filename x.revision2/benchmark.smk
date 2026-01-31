rule all:
    input:
        expand("benchmark/{tool}/COLO829_{purity}.60x.COLO829_BL.60x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.60x.COLO829_BL.45x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.60x.COLO829_BL.30x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.45x.COLO829_BL.45x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.45x.COLO829_BL.30x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.45x.COLO829_BL.15x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.30x.COLO829_BL.30x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.30x.COLO829_BL.15x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),
        expand("benchmark/{tool}/COLO829_{purity}.15x.COLO829_BL.15x/summary.txt", tool = ["deepsomatic", "clairS"], purity = ["10","20","40","60","80","100"]),

ClairS_sif = "/mnt/backedup/home/jiaZ/working/containers/sif/clairs/clairs_v0.1.7.sif"
rule benchmark_snvs:
    input:
        snv = "simulated/snvs/{tool}/COLO829_{purity}.{depth_t}.COLO829_BL.{depth_n}/output.passed.vcf.gz",
        genome = "genome/chr22.fasta"
    output:
        "benchmark/{tool}/COLO829_{purity}.{depth_t}.COLO829_BL.{depth_n}/summary.txt"
    params:
        outdir = "benchmark/{tool}/COLO829_{purity}.{depth_t}.COLO829_BL.{depth_n}",
        gs = "genome/chr22.vcf.gz",
    threads: 8
    resources:
        mem = 24,
        walltime = 2
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        echo -e "Filter\tType\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output}
        
        echo -ne "PASS\t" >> {output}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {input.snv} \
        --input_filter_tag 'PASS' \
        --threads {threads} \
        --ref_fn {input.genome} \
        --output_dir {params.outdir} --roc_fn {params.outdir}/roc | grep -E 'SNV' >> {output}
        """