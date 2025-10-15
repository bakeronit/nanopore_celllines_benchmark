configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"

rule all:
    input:
        #expand("analysis/bam/{sample}.bam", sample=['COLO829','COLO829_BL', "HCC1937", "HCC1937_BL"])
        expand("analysis/snvs/clairS/{sample_t}.{sample_n}/output.vcf.gz", zip, sample_t=["COLO829","HCC1937"], sample_n=["COLO829_BL","HCC1937_BL"]),
        expand("analysis/snvs/deepsomatic/COLO829.COLO829_BL/output.{chrom}.vcf.gz", chrom=[f'chr{i}' for i in range(1,23)] + ['chrX', 'chrY', 'chrM']),
        expand("analysis/snvs/deepsomatic/HCC1937.HCC1937_BL/output.{chrom}.vcf.gz", chrom=[f'chr{i}' for i in range(1,23)] + ['chrX', 'chrY', 'chrM']),
        expand("analysis/svs/delly/{sample_t}.{sample_n}/{sample_t}.{sample_n}.vcf", zip, sample_t=["COLO829","HCC1937"], sample_n=["COLO829_BL","HCC1937_BL"]),
        expand("analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.pass.svtype.txt", zip, sample_t=["COLO829","HCC1937"], sample_n=["COLO829_BL","HCC1937_BL"]),
        #expand("analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf", zip, sample_t=["COLO829","HCC1937"], sample_n=["COLO829_BL","HCC1937_BL"]),
        expand("analysis/svs/severus/{sample_t}.{sample_n}/somatic_SVs/severus_somatic.vcf", zip, sample_t=["COLO829","HCC1937"], sample_n=["COLO829_BL","HCC1937_BL"])

wildcard_constraints:
    sample = "|".join(['COLO829','COLO829_BL', 'HCC1937', 'HCC1937_BL']),
    sample_t = "|".join(['COLO829','HCC1937']),
    sample_n = "|".join(['COLO829_BL','HCC1937_BL']),

rule align_t2t_minimap2:
    input:
        bam = "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/bam/R10/sup/{sample}.bam",
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa"
    output:
        bam = "analysis/bam/{sample}.bam",
        bai = "analysis/bam/{sample}.bam.bai"
    envmodules:
        "samtools/1.17",
        "minimap2/2.26"
    threads: 24
    resources:
        mem = 36,
        walltime = 48
    shell:
        """
        samtools fastq -@8 -T"*" {input.bam} | \
        minimap2 -y --MD -ax map-ont -t 8 {input.genome} - | \
        samtools sort -@8 -O BAM --write-index -o {output.bam}##idx##{output.bai}
        """

def estimated_mem(chrom):
    chrom_sizes = {
    'chr1': 248387328,
    'chr2': 242696752,
    'chr3': 201105948,
    'chr4': 193574945,
    'chr5': 182045439,
    'chr6': 172126628,
    'chr7': 160567428,
    'chr8': 146259331,
    'chr9': 150617247,
    'chr10': 134758134,
    'chr11': 135127769,
    'chr12': 133324548,
    'chr13': 113566686,
    'chr14': 101161492,
    'chr15': 99753195,
    'chr16': 96330374,
    'chr17': 84276897,
    'chr18': 80542538,
    'chr19': 61707364,
    'chr20': 66210255,
    'chr21': 45090682,
    'chr22': 51324926,
    'chrM': 16569,
    'chrX': 154259566,
    'chrY': 62460029
}
    min_mem = 30
    base_chr = "chr22"
    base_mem = 70
    base_size = chrom_sizes[base_chr]
    if not chrom_sizes.get(chrom):
        return min_mem
    mem = int(chrom_sizes[chrom] * base_mem * 1.2 / base_size) # allow 20% buffer
    if mem > min_mem:
        return mem
    return min_mem

## call somatic snv and indels
include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/snv_calling/clairS.smk"
use rule call_somatic_snv_clairS as call_somatic_snv_clairS_t2t with:
    input:
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa",
        tumor_bam = "analysis/bam/{sample_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai"
    output:
        "analysis/snvs/clairS/{sample_t}.{sample_n}/output.vcf.gz"
    params:
        platform = "ont_r10_dorado_sup_5khz",
        outdir = "analysis/snvs/clairS/{sample_t}.{sample_n}",
        clair3_model = "/mnt/backedup/home/jiaZ/working/data/ont_models/clair3_models/r1041_e82_400bps_sup_v430",
        indel_option = "--enable_indel_calling"
    log:
        "logs/clairS/{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/clairS/{sample_t}.{sample_n}.benchmark.txt"

include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/snv_calling/deepsomatic.smk" # just import the sif file path
rule call_somatic_snv_deepsomatic_t2t:
    input:
        genome="/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa",
        tumor_bam = "analysis/bam/{sample_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.{chrom}.vcf.gz"
    params:
        model="/mnt/backedup/home/jiaZ/working/data/ont_models/dpsomatic_model/weights-143-0.987994.ckpt"
    threads: 24
    envmodules:
        "singularity/3.7.1"
    resources:
        mem=lambda wildcards: estimated_mem(wildcards.chrom),
        walltime=40
    shell:
        """
        singularity exec {DP_somatic_sif} /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --ref={input.genome} \
        --model_type=ONT_R104 \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --num_shards={threads} \
        --logging_dir=logs/deepsomatic/{wildcards.sample_t}.{wildcards.sample_n} \
        --customized_model={params.model} \
        --regions={wildcards.chrom}
        """

#### call somatic SV
delly = config['delly']['path']
rule call_somatic_sv_delly_t2t:
    input:
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa",
        tumour_bam = "analysis/bam/{sample_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai"
    output:
        bcf = "analysis/svs/delly/{sample_t}.{sample_n}/{sample_t}.{sample_n}.pre.bcf",
        samples_tsv = "analysis/svs/delly/{sample_t}.{sample_n}/{sample_t}.{sample_n}.tsv",
        filtered_vcf = "analysis/svs/delly/{sample_t}.{sample_n}/{sample_t}.{sample_n}.vcf"
    params:
        excl = f"-x {config['delly']['bed']}" if config['delly']['bed'] != None else ""
    threads: 2
    envmodules:
        "samtools/1.17",
        "bcftools/1.16"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        export OMP_NUM_THREADS={threads} # delly primarily parallises on the sample level.
        {delly} lr \
            -t ALL \
            -y ont \
            -o {output.bcf} {params.excl} \
            -g {input.genome} {input.tumour_bam} {input.normal_bam}
        
        printf "{wildcards.sample_t}\\ttumor\\n{wildcards.sample_n}\\tcontrol\\n" > {output.samples_tsv}

        {delly} filter \
            -f somatic \
            -p \
            -s {output.samples_tsv} \
            {output.bcf} | bcftools view -Ov > {output.filtered_vcf}
        """

include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/sv_calling/nanomonsv.smk"
use rule nanomonsv_parse as nanomonsv_parse_t2t with:
    input:
        bam = "analysis/bam/{sample}.bam",
        bai = "analysis/bam/{sample}.bam.bai",
    output:
        multiext("analysis/svs/nanomonsv/{sample}/{sample}.",
                "bp_info.sorted.bed.gz", 
                "bp_info.sorted.bed.gz.tbi",
                "deletion.sorted.bed.gz",
                "insertion.sorted.bed.gz",
                "rearrangement.sorted.bedpe.gz")
    params:
        prefix="analysis/svs/nanomonsv/{sample}/{sample}"

rule call_somatic_sv_nanomonsv_get_t2t:
    input:
        multiext("analysis/svs/nanomonsv/{sample_t}/{sample_t}.","bp_info.sorted.bed.gz", "bp_info.sorted.bed.gz.tbi","deletion.sorted.bed.gz","insertion.sorted.bed.gz","rearrangement.sorted.bedpe.gz"),
        multiext("analysis/svs/nanomonsv/{sample_n}/{sample_n}.","bp_info.sorted.bed.gz","bp_info.sorted.bed.gz.tbi","deletion.sorted.bed.gz","insertion.sorted.bed.gz","rearrangement.sorted.bedpe.gz"),        
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa",
        tumour_bam = "analysis/bam/{sample_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai",
        control_panel_path = config['nanomonsv']['control_panel_path'],
    output:
        txt = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.txt",
        vcf = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.vcf",
        sbnd_txt = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd.result.txt",
        sread = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.supporting_read.txt"
    params:
        tumour_prefix = "analysis/svs/nanomonsv/{sample_t}/{sample_t}",
        normal_prefix = "analysis/svs/nanomonsv/{sample_n}/{sample_n}",
        panel_prefix="hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control",
        final_outdir = "analysis/svs/nanomonsv/{sample_t}.{sample_n}"
    threads: 20
    resources:
        mem = 48,
        walltime = 48
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        nanomonsv get {params.tumour_prefix} {input.tumour_bam} {input.genome} \
        --control_prefix {params.normal_prefix} --control_bam {input.normal_bam} \
        --single_bnd --use_racon --min_indel_size 10 --qv20 \
        --control_panel_prefix {input.control_panel_path}/{params.panel_prefix} --processes {threads}

        mkdir -p {params.final_outdir}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}/{wildcards.sample_t}.nanomonsv.result.txt {output.txt}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}/{wildcards.sample_t}.nanomonsv.result.vcf {output.vcf}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}/{wildcards.sample_t}.nanomonsv.sbnd.result.txt {output.sbnd_txt}
        mv analysis/svs/nanomonsv/{wildcards.sample_t}/{wildcards.sample_t}.nanomonsv.supporting_read.txt {output.sread}
        """

use rule nanomonsv_filter_simple_repeat_svtype as nanomonsv_filter_simple_repeat_svtype_t2t with:
    input:
        result = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.txt",
        simple_repeat = config['nanomonsv']['simple_repeat']
    output:
        filt = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.txt",
        passed = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.pass.txt",
        svtype = "analysis/svs/nanomonsv/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.pass.svtype.txt"

include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/sv_calling/savana.smk"
use rule call_somatic_sv_savana as call_somatic_sv_savana_t2t with:
    input:
        tumour_bam = "analysis/bam/{sample_t}.bam",
        tumour_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai",
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa"
    output:
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.bedpe",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints_read_support.tsv",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.sv_breakpoints.vcf",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf"
    params:
        outdir = "analysis/svs/savana/{sample_t}.{sample_n}"
    log: 
        "logs/savana/{sample_t}.{sample_n}.log"
    benchmark: 
        "benchmarks/savana/{sample_t}.{sample_n}.benchmark.txt"


include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/snv_calling/phasing.smk"
include: "../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/sv_calling/severus.smk"
Clair3_sif = config['clair3']['sif']
rule call_germline_snv_clair3_t2t:
    input:
        bam = "analysis/bam/{sample}.bam",
        bai = "analysis/bam/{sample}.bam.bai",
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa"
    output:
        "analysis/snvs/clair3/{sample}/phased_merge_output.vcf.gz",
    params:
        outdir = "analysis/snvs/clair3/{sample}",
        model = "/mnt/backedup/home/jiaZ/working/data/ont_models/clair3_models/r1041_e82_400bps_sup_v430"
    threads: 24
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {Clair3_sif} /opt/bin/run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.genome} \
            --sample_name={wildcards.sample} \
            --threads={threads} \
            --platform="ont" \
            --model_path={params.model} \
            --enable_phasing \
            --include_all_ctgs \
            --remove_intermediate_dir \
            --output={params.outdir}
        """

use rule haplotagging_whatshap as haplotagging_whatshap_t2t with:
    input:
        genome = "/mnt/backedup/home/jiaZ/working/data/genome/chm13/chm13v2.0.fa",
        vcf = "analysis/snvs/clair3/{sample}/phased_merge_output.vcf.gz",
        bam = "analysis/bam/{sample}.bam",
    output:
        bam = "analysis/bam/{sample}.haplotagged.bam"
    params:
        outdir = "analysis/bam/{sample}"
    log:
        "logs/whatshap/{sample}.log"
    benchmark:
        "benchmarks/whatshap/{sample}.benchmark.txt"

use rule index_haplotagged_bam as index_haplotagged_bam_t2t with:
    input:
        "analysis/bam/{sample}.haplotagged.bam"
    output:
        "analysis/bam/{sample}.haplotagged.bam.bai"

use rule call_somatic_sv_severus as call_somatic_sv_severus_t2t with:
    input:
        hp_tagged_tumour_bam = "analysis/bam/{sample_t}.haplotagged.bam",
        hp_tagged_normal_bam = "analysis/bam/{sample_n}.haplotagged.bam",
        hp_tagged_tumour_bai = "analysis/bam/{sample_t}.haplotagged.bam.bai",
        hp_tagged_normal_bai = "analysis/bam/{sample_n}.haplotagged.bam.bai",
        phased_vcf = "analysis/snvs/clair3/{sample_n}/phased_merge_output.vcf.gz",
        vntr_bed = config['severus']['vntr']
    output:
        "analysis/svs/severus/{sample_t}.{sample_n}/somatic_SVs/severus_somatic.vcf"
    params:
        outdir = "analysis/svs/severus/{sample_t}.{sample_n}"
    log:
        "logs/severus/{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/severus/{sample_t}.{sample_n}.txt"
    resources:
        mem = 80,
        walltime = 48
