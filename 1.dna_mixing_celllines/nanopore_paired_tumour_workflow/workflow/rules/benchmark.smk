rule benchmark_somatic_clairs:
    input:
        snv = "analysis/snvs/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/output.vcf.gz",
        ref = config['reference']['file'] 
    output:
        "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/fn.vcf",
        "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/fp_fn.vcf",
        "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/fp.vcf",
        "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/tp.vcf",
        txt = "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}/summary.txt"
    params:
        output_dir = "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}",
        output_pass_dir = "analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{sample_t}.{sample_n}_pass",
        gs = lambda w: config['gold_standard']['somatic'][w.sample_t]
    threads: 4
    resources:
        mem = 10,
        walltime = 1
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        echo -e "Filter\tType\tPrecision\tRecall\tF1-score\tTP\tFP\tFN" > {output.txt}
        
        echo -ne "ALL\t" >> {output.txt}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {input.snv} \
        --threads {threads} \
        --ref_fn {input.ref} \
        --output_dir {params.output_dir} | grep -E 'SNV' >> {output.txt}
        
        echo -ne "PASS\t" >> {output.txt}
        singularity exec {ClairS_sif} python /opt/bin/clairs.py compare_vcf \
        --truth_vcf_fn {params.gs} \
        --input_vcf_fn {input.snv} \
        --input_filter_tag 'PASS' \
        --threads {threads} \
        --ref_fn {input.ref} \
        --output_dir {params.output_pass_dir} | grep -E 'SNV' >> {output.txt}
        """

use rule benchmark_somatic_clairs as benchmark_germline_clair3 with:
    input:
        snv = "analysis/snvs/{caller}/{flowcell}/{mode}/{sample}/merge_output.vcf.gz",
        ref = config['reference']['file']
    output:
        "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}/fn.vcf",
        "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}/fp_fn.vcf",
        "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}/fp.vcf",
        "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}/tp.vcf",
        txt = "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}/summary.txt"
    params:
        gs = lambda w: config['gold_standard']['germline'][w.sample],
        output_dir = "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}",
        output_pass_dir = "analysis/benchmarks/snvs/germline/{caller}/{flowcell}/{mode}/{sample}_pass"

use rule benchmark_germline_clair3 as benchmarks_germline_pp_dpvariant with:
    input:
        snv = "analysis/snvs/{caller}/{flowcell}/{mode}/{sample}/{sample}.vcf.gz",
        ref = config['reference']['file']
