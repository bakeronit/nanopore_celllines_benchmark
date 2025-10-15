from pathlib import Path

snv_path=Path("/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/x.revision/analysis/snvs")


rule all:
    input:
        expand(snv_path / "{tool}/COLO829.COLO829_BL/benchmark/summary.txt", tool=['clairS', 'deepsomatic']),
        expand(snv_path / "{tool}/HCC1937.HCC1937_BL/benchmark/summary.txt", tool=['clairS', 'deepsomatic']),
        expand(snv_path / "{tool}/COLO829.COLO829_BL/hg38/output.snv.hg38_to_t2t.vcf.gz", tool=['clairS', 'deepsomatic']),
        expand(snv_path / "{tool}/HCC1937.HCC1937_BL/hg38/output.snv.hg38_to_t2t.vcf.gz", tool=['clairS', 'deepsomatic']),


wildcard_constraints:
    sample_t="COLO829|HCC1937",
    sample_n="COLO829_BL|HCC1937_BL",
    tool="clairS|deepsomatic"


rule get_passed_snv_clairs:
    input:
        snv_path / "clairS/{sample_t}.{sample_n}/output.vcf.gz"
    output:
        snv_path / "clairS/{sample_t}.{sample_n}/output.snv.vcf.gz",
    envmodules:
        "bcftools/1.19",
        "htslib/1.22.1"
    threads: 1
    resources:
        mem=1,
        walltime=1
    shell:
        """
        bcftools view -f PASS {input} |bgzip > {output}
        tabix -p vcf {output}
        """        

rule liftover_t2t:
    input:
        snv = snv_path / "{tool}/{sample_t}.{sample_n}/output.snv.vcf.gz",
        reference = "/mnt/backedup/home/jiaZ/working/data/genome/reference.fasta"
    output:
        snv = snv_path / "{tool}/{sample_t}.{sample_n}/output.snv.t2t_to_hg38.vcf.gz",
        snv_reject = snv_path / "{tool}/{sample_t}.{sample_n}/output.snv.t2t_to_hg38.rejected.vcf.gz",
    params:
        chain = "/mnt/lustre/working/lab_nicw/jiaZ/bioprojects/nanopore_celllines_benchmark/x.revision/data/chm13v2-grch38.chain"
    envmodules:
        "gatk/4.6.2.0"
    threads: 2
    resources:
        mem=24,
        walltime=1
    shell:
        """
        gatk --java-options "-Xms4g -Xmx24g -Djava.io.tmpdir=/scratch" \
        LiftoverVcf \
        --CHAIN {params.chain} \
        --INPUT {input.snv} \
        --OUTPUT {output.snv} \
        --REFERENCE_SEQUENCE {input.reference} \
        --REJECT {output.snv_reject}
        """

ClairS_sif="/mnt/backedup/home/jiaZ/working/imgs/clairs/clairs_v0.1.7.sif"
rule benchmark_snv:
    input:
        snv = snv_path / "{tool}/{sample_t}.{sample_n}/output.snv.t2t_to_hg38.vcf.gz",
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/reference.fasta"
    output:
        snv_path / "{tool}/{sample_t}.{sample_n}/benchmark/summary.txt"
    params:
        outdir = str(snv_path) + "/{tool}/{sample_t}.{sample_n}/benchmark",
        gs = lambda wildcards: f"/mnt/backedup/home/jiaZ/working/general/goldstandard/vcfs/{wildcards.sample_t.lower()}/merged_normed_isec_snv.hom100.goldstandard.vcf.gz"
    threads: 8
    resources:
        mem = 20,
        walltime = 1
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
        --output_dir {params.outdir} | grep -E 'SNV' >> {output}
        """


hg38_snv_path = Path("/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/snvs")
rule liftover_hg38:
    input:
        snv = hg38_snv_path / "{tool}/R10/sup/{sample_t}.{sample_n}/output.vcf.gz",
        reference = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa"
    output:
        filtered_hg38_snv = snv_path / "{tool}/{sample_t}.{sample_n}/hg38/output.snv.vcf.gz",
        snv = snv_path / "{tool}/{sample_t}.{sample_n}/hg38/output.snv.hg38_to_t2t.vcf.gz",
        snv_reject = snv_path / "{tool}/{sample_t}.{sample_n}/hg38/output.snv.hg38_to_t2t.rejected.vcf.gz",
    params:
        chain = "/mnt/lustre/working/lab_nicw/jiaZ/bioprojects/nanopore_celllines_benchmark/x.revision/data/grch38-chm13v2.chain"
    envmodules:
        "gatk/4.6.2.0",
        "bcftools/1.19",
        "htslib/1.22.1"
    threads: 2
    resources:
        mem=24,
        walltime=1
    shell:
        """
        bcftools view -v snps -f PASS {input.snv} | bgzip > {output.filtered_hg38_snv}
        tabix -p vcf {output.filtered_hg38_snv}
        
        gatk --java-options "-Xms4g -Xmx24g -Djava.io.tmpdir=/scratch" \
        LiftoverVcf \
        --CHAIN {params.chain} \
        --INPUT {output.filtered_hg38_snv} \
        --OUTPUT {output.snv} \
        --REFERENCE_SEQUENCE {input.reference} \
        --REJECT {output.snv_reject}
        """
        