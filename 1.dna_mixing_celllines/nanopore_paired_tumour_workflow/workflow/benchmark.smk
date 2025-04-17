workdir: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work"

##### setup report #####
#configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/config.hcc1937.yaml"
configfile: "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config/bc.config.yaml"

include: "rules/common.smk"
include: "rules/snv_calling/clair3.smk"
include: "rules/snv_calling/clairS.smk"
include: "rules/snv_calling/deepsomatic.smk"
include: "rules/snv_calling/deepvariant.smk"

pairs = generate_paired_samples(samples_df)
[f"analysis/benchmark/snvs/clairS/somatic/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE]

rule all:
    input:
        #"analysis/benchmark/snvs/germline/R10/sup/clair3/COLO829_BL/summary.txt",
        #"analysis/benchmark/snvs/germline/R10/sup/deepvariant/COLO829_BL/summary.txt",
        #"analysis/benchmark/snvs/germline/R10/sup/clair3/HCC1937_BL/summary.txt",
        #"analysis/benchmark/snvs/germline/R10/sup/deepvariant/HCC1937_BL/summary.txt",
        #[f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/clairS/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE],
        #[f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/deepsomatic/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs for m in MODE]
        ## benchmark without chrXY, just for HCC1937
        [f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/clairS/{pair['tumour']}.{pair['normal']}_noXY/summary.txt" for pair in pairs for m in MODE],
        [f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/deepsomatic/{pair['tumour']}.{pair['normal']}_noXY/summary.txt" for pair in pairs for m in MODE],
        #[f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/clairS/{pair['tumour']}.{pair['normal']}_noXY_passed/roc" for pair in pairs for m in MODE],
        #[f"analysis/benchmark/snvs/somatic/{pair['flowcell_version']}/{m}/deepsomatic/{pair['tumour']}.{pair['normal']}_noXY_passed/roc" for pair in pairs for m in MODE]


def get_snv_results(wildcards, exclude_xy=False):
    match wildcards.tool:
        case "clair3":
            return rules.call_germline_snv_clair3.output
        case "deepvariant":
            return rules.call_germline_snv_deepvariant.output.vcf
        case "clairS":
            if exclude_xy:
                return "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/output.no_chrXY.vcf.gz"
            return rules.call_somatic_snv_clairS.output
        case "deepsomatic":
            if exclude_xy:
                return "analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.somatic.no_chrXY.vcf.gz"
            return rules.get_somatic_deepsomatic.output.vcf
        case _:
            raise VauleError("tool can only be clair3, deepvariant for germline, and clairS, deepsomatic for somatic!")

def get_snv_goldstandard(sample_id):
    cellline = samples_df[samples_df['sample_id'] == sample_id].donor_id.unique().tolist()[0]
    type = 'somatic' if samples_df[samples_df['sample_id'] == sample_id].type.unique().tolist()[0] == "tumour" else "germline"
    return config['snv'][type][cellline]

rule get_somatic_deepsomatic:
    input:
        rules.call_somatic_snv_deepsomatic.output
    output:
        vcf="analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.somatic.vcf.gz",
        tbi="analysis/snvs/deepsomatic/{flowcell}/{mode}/{sample_t}.{sample_n}/output.somatic.vcf.gz.tbi"
    threads: 1
    resources:
        mem = 2,
        walltime = 1
    envmodules:
        "bcftools/1.16",
        "htslib/1.17"
    shell:
        """
        bcftools filter -i 'GT="1/1"' {input} |bgzip > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule benchmark_germline_snv_calling:
    input:
        snv = get_snv_results,
        genome = config['reference']['file'],
    output:
        "analysis/benchmark/snvs/germline/{flowcell}/{mode}/{tool}/{sample}/summary.txt"
    params:
        outdir = "analysis/benchmark/snvs/germline/{flowcell}/{mode}/{tool}/{sample}",
        passed_outdir = "analysis/benchmark/snvs/germline/{flowcell}/{mode}/{tool}/{sample}_passed",
        gs = lambda w: get_snv_goldstandard(w.sample)
    threads: 8
    resources:
        mem = 50,
        walltime = 2
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        mkdir -p {params.passed_outdir}
        echo -e "Filter\tType\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output}
        
        echo -ne "ALL\t" >> {output}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {input.snv} \
        --threads {threads} \
        --ref_fn {input.genome} \
        --output_dir {params.outdir} | grep -E 'SNV' >> {output}
        
        echo -ne "PASS\t" >> {output}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {input.snv} \
        --input_filter_tag 'PASS' \
        --threads {threads} \
        --ref_fn {input.genome} \
        --output_dir {params.passed_outdir} --roc_fn {params.passed_outdir}/roc | grep -E 'SNV' >> {output}
        """

use rule benchmark_germline_snv_calling as benchmark_somatic_snv_calling with:
    output:
        "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}/summary.txt"
    params:
        outdir = "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}",
        passed_outdir = "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}_passed",
        gs = lambda w: get_snv_goldstandard(w.sample_t)


# benchmark snv calling with out chrXY for HCC1937
use rule benchmark_germline_snv_calling as benchmark_somatic_snv_calling_noXY with:
    input:
        snv = lambda w: get_snv_results(w, exclude_xy=True),
        genome = config['reference']['file'],
    output:
        "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}_noXY/summary.txt",
    params:
        outdir = "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}_noXY",
        passed_outdir = "analysis/benchmark/snvs/somatic/{flowcell}/{mode}/{tool}/{sample_t}.{sample_n}_noXY_passed",
        gs = "/mnt/backedup/home/jiaZ/working/general/goldstandard/vcfs/hcc1937/merged_normed_isec_snv.hom100.goldstandard.noXY.vcf.gz"
