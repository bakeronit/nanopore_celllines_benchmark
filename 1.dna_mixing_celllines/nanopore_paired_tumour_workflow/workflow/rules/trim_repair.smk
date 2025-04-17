# trim the heading of 40bp and tailing of 20bp for each read.
# after trimmming, need to fix the MM tag.

modkit = config['modkit']
rule trim_reads:
    input:
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam"
    output:
        temp("analysis/ubam/{flowcell}/{mode}/{sample}/{run}.trimmed.sorted.ubam")
    threads: 12
    resources:
        mem = 20,
        walltime = 24
    envmodules:
        "samtools/1.17"
    params:
        headcrop = config['filter']['headcrop'],
        tailcrop = config['filter']['tailcrop'],
        minlen = config['filter']['minlen'],
    shell:
        """
        samtools fastq -T"MM,ML,qs" {input} | \
        chopper --headcrop {params.headcrop} --tailcrop {params.tailcrop} --minlength {params.minlen} --threads {threads} | \
        samtools import -T"*" - | \
        samtools sort -@{threads} -n > {output}
        """

rule sorted_basecalling_ubam:
    input:
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam"
    output:
        temp("analysis/ubam/{flowcell}/{mode}/{sample}/{run}.sorted.ubam")
    threads: 12
    resources:
        mem = 20,
        walltime = 10
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools sort -@{threads} -n {input} > {output}
        """

rule repair_MMtag:
    input:
        original = "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.sorted.ubam",
        trimmed = "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.trimmed.sorted.ubam"
    output:
        temp("analysis/ubam/{flowcell}/{mode}/{sample}/{run}.trimmed_repaired.ubam")
    resources:
        mem = 20,
        walltime = 10
    threads: 1
    shell:
        """
        {modkit} repair -d {input.original} -a {input.trimmed} -o {output}
        """
