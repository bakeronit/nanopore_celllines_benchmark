#samples = [f"COLO829_{purity}" for purity in range(40,100,20)] + ["COLO829","COLO829_BL"]
samples = [f"HCC1937_{purity}" for purity in range(40,100,20)] + ["HCC1937","HCC1937_BL"]
wildcard_constraints:
    sample = "|".join(samples),
    depth = "|".join(["15x","30x","45x","60x"])

rule all:
    input:
        expand("analysis/bam/{sample}.{depth}.mosdepth.summary.txt",sample=samples[:-1], depth=["15x","30x","45x","60x"]),
        #expand("analysis/bam/COLO829_BL.{depth}.mosdepth.summary.txt",depth=["15x","30x","45x","60x"])
        expand("analysis/bam/HCC1937_BL.{depth}.mosdepth.summary.txt",depth=["15x","30x","45x","60x"])

rule subsample_bam:
    input:
        mosdepth = "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/{sample}.mosdepth.summary.txt",
        bam = "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/bam/R10/sup/{sample}.bam",
        bai = "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/bam/R10/sup/{sample}.bam.bai"
    output:
        bam = "analysis/bam/{sample}.{depth}.bam",
        bai = "analysis/bam/{sample}.{depth}.bam.bai"
    params:
        t_depth = lambda w: int(w.depth.replace('x',''))  # get the integer from wildcards.
    threads: 10
    resources:
        mem = 10,
        walltime = 8
    log:
        "logs/subsample/{sample}.{depth}.log"
    envmodules:
        "samtools/1.17"
    shell:
        """
        ave_depth=`grep total {input.mosdepth} | cut -f4`  # the genome-wide average sequencing depth
        if (( $(bc <<< "scale=2; {params.t_depth}/$ave_depth >= 1 ") )); then  # there are some sample didn't reach 60x but close to, so just use the original BAM.
            echo "{wildcards.sample} sequencing depth $ave_depth didn't reach the target depth {params.t_depth}, using symlink instead." > {log}
            ln -s {input.bam} {output.bam}
            ln -s {input.bai} {output.bai}
        else
            prop=`echo "scale=2; {params.t_depth}/$ave_depth" | bc`
            echo "subsampling $prop from {wildcards.sample} with $ave_depth to reach {params.t_depth}" > {log}
            samtools view -b -@{threads} --subsample $prop --subsample-seed 2024 --write-index -o {output.bam}##idx##{output.bai} {input.bam}
        fi
        """

rule check_depth:
    input:
        "analysis/bam/{sample}.{depth}.bam",
    output:
        "analysis/bam/{sample}.{depth}.mosdepth.global.dist.txt",
        "analysis/bam/{sample}.{depth}.mosdepth.summary.txt",
    params:
        prefix = "analysis/bam/{sample}.{depth}"
    threads: 10
    resources:
        mem = 10,
        walltime = 8
    envmodules:
        "mosdepth/0.2.9"
    shell:
        """
        mosdepth -n -t {threads} {params.prefix} {input}
        """

