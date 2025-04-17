workdir: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work"

##### setup report #####
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/bc.config.yaml"

include: "rules/common.smk"
include: "rules/snv_calling/clair3.smk"
include: "rules/snv_calling/clairS.smk"
include: "rules/snv_calling/deepsomatic.smk"
include: "rules/snv_calling/deepvariant.smk"

pairs = generate_paired_samples(samples_df)
[f"analysis/benchmark/indels/clairS/somatic/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE]

rule all:
    input:
        #"analysis/benchmark/indels/germline/R10/sup/clair3/COLO829_BL/summary.txt",
        #"analysis/benchmark/indels/germline/R10/sup/deepvariant/COLO829_BL/summary.txt",
        #"analysis/benchmark/indels/germline/R10/sup/clair3/HCC1937_BL/summary.txt",
        #"analysis/benchmark/indels/germline/R10/sup/deepvariant/HCC1937_BL/summary.txt",
        [f"analysis/benchmark/indels/somatic/{pair['flowcell_version']}/{m}/clairS/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE],
        [f"analysis/benchmark/indels/somatic/{pair['flowcell_version']}/{m}/deepsomatic/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE]


def get_vcf_results(wildcards):
    match wildcards.tool:
        case "clair3":
            return rules.call_germline_snv_clair3.output
        case "deepvariant":
            return rules.call_germline_snv_deepvariant.output.vcf
        case "clairS":
            return "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/indel.vcf.gz"
        case "deepsomatic":
            return "analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.somatic.vcf.gz"
        case _:
            raise VauleError("tool can only be clair3, deepvariant for germline, and clairS, deepsomatic for somatic!")

def get_indel_goldstandard(sample_id):
    cellline = samples_df[samples_df['sample_id'] == sample_id].donor_id.unique().tolist()[0]
    type = 'somatic' if samples_df[samples_df['sample_id'] == sample_id].type.unique().tolist()[0] == "tumour" else "germline"
    return config['indel'][type][cellline]

rule benchmark_germline_indel_calling:
    input:
        vcf = get_vcf_results,
        genome = config['reference']['file'],
    output:
        "analysis/benchmark/indels/germline/{flowcell}/{mode}/{tool}/{sample}/summary.txt"
    params:
        outdir = "analysis/benchmark/indels/germline/{flowcell}/{mode}/{tool}/{sample}",
        norm_vcf = "analysis/snvs/{tool}/{flowcell}/{mode}/{sample}/norm_indel.vcf.gz",
        gs = lambda w: get_indel_goldstandard(w.sample)
    threads: 2
    resources:
        mem = 10,
        walltime = 2
    envmodules:
        "singularity/3.7.1",
        "bcftools/1.19",
        "htslib/1.19.1"
    shell:
        """
        bcftools view -v indels {input.vcf} | bcftools norm -a --atom-overlaps . | bgzip > {params.norm_vcf} 
        echo -e "Type\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output}
        
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {params.norm_vcf} \
        --input_filter_tag 'PASS' \
        --benchmark_indel \
        --threads {threads} \
        --ref_fn {input.genome} \
        --output_dir {params.outdir} | grep -E 'INDEL|INS|DEL' | awk '{{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7}}' >> {output}
        """

use rule benchmark_germline_indel_calling as benchmark_somatic_indel_calling with:
    output:
        "analysis/benchmark/indels/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}/summary.txt"
    params:
        outdir = "analysis/benchmark/indels/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}",
        norm_vcf = "analysis/snvs/{tool}/{flowcell}/{mode}/{sample_t}.{sample_n}/norm_indel.vcf.gz",
        gs = lambda w: get_indel_goldstandard(w.sample_t)

