workdir: "gs"
configfile: "/path/to/config.yaml"


wildcard_constraints:
    cell = "|".join(['colo829','hcc1937']),
    tumor = "|".join(['colo829','hcc1937']),
    normal = "|".join(['colo829_bl','hcc1937_bl']),


rule all:
    input:
        #expand("analysis/svs/jasmine_merge/colo829/{lib}/colo829.colo829_bl.{tool}.simple.vcf", lib=['B','C','D'], tool=['gridss','lumpy','delly']),
        #expand("analysis/svs/jasmine_merge/hcc1937/{lib}/hcc1937.hcc1937_bl.{tool}.simple.vcf", lib=['A','B','C'], tool=['gridss','lumpy','delly']),
        #expand("analysis/svs/jasmine_merge/colo829/{lib}/colo829.colo829_bl.merged.vcf",lib=['B','C','D']),
        #expand("analysis/svs/jasmine_merge/hcc1937/{lib}/hcc1937.hcc1937_bl.merged.vcf",lib=['A','B','C']),
        expand("analysis/svs/jasmine_merge/hcc1937/hcc1937.final_merged.vcf"),
        expand("analysis/svs/jasmine_merge/colo829/colo829.final_merged.vcf")

rule get_simple_type:
    input:
        gridss = "analysis/svs/gridss/{cell}/{lib}/{tumor}.{normal}/{tumor}.bam.gripss.filtered.vcf.gz",
        lumpy = "analysis/svs/lumpy/{cell}/{lib}/{tumor}.{normal}.somatic.vcf",
        delly = "analysis/svs/delly/{cell}/{lib}/{tumor}.{normal}.vcf",
    output:
        gridss = "analysis/svs/jasmine_merge/{cell}/{lib}/{tumor}.{normal}.gridss.simple.vcf",
        lumpy = "analysis/svs/jasmine_merge/{cell}/{lib}/{tumor}.{normal}.lumpy.simple.vcf",
        delly = "analysis/svs/jasmine_merge/{cell}/{lib}/{tumor}.{normal}.delly.simple.vcf",
    params:
        jasmine_env = config['jasmine_env'],
        script = config['simple_script'],
    shell:
        """
        module load conda-envs/base bcftools/1.19
        conda activate {params.jasmine_env}
        
        python {params.script} {input.gridss} -t gridss | bcftools view -f 'PASS,.' | bcftools filter -i '( (SVTYPE="DUP" || SVTYPE="DEL" || SVTYPE="INV") && SVLEN>=50 ) || (SVTYPE="INS") || (SVTYPE="TRA")' > {output.gridss}
        python {params.script} {input.delly} -t delly | bcftools view -f 'PASS,.' | bcftools filter -i '( (SVTYPE="DUP" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="INS") && SVLEN>=50 ) || (SVTYPE="TRA")' > {output.delly}
        python {params.script} {input.lumpy} -t lumpy | bcftools view -f 'PASS,.' | bcftools filter -i '( (SVTYPE="DUP" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="INS") && SVLEN>=50 ) || (SVTYPE="TRA")' > {output.lumpy}    
        """

rule jasmine_merge_tools:
    """
    Merge from 3 tools and get SVs present in at least 2 tools. For better merging, the DUP was converted into INS with noting in INFO with OLDTYPE=DUP. \
    I didn't convert them back since I also need to merge SVs events called in 3 libraries.
    """
    input:
        gridss = rules.get_simple_type.output.gridss,
        lumpy = rules.get_simple_type.output.lumpy,
        delly = rules.get_simple_type.output.delly,
        genome = config['reference']
    output:
        merged = "analysis/svs/jasmine_merge/{cell}/{lib}/{tumor}.{normal}.merged.vcf",
        overlap = "analysis/svs/jasmine_merge/{cell}/{lib}/{tumor}.{normal}.merged.supp2.vcf"
    params:
        jasmine_env = config['jasmine_env'],
        out_dir = "analysis/svs/jasmine_merge/{cell}/{lib}/jasmine_out",
        simple_filelist = "analysis/svs/jasmine_merge/{cell}/{lib}.simple.files.txt"
    shell:
        """
        module load conda-envs/base
        module load bcftools/1.19
        conda activate {params.jasmine_env}

        ls {input.delly} {input.gridss} {input.lumpy} > {params.simple_filelist}
        jasmine --preprocess_only --pre_normalize --dup_to_ins file_list={params.simple_filelist} out_dir={params.out_dir} genome_file={input.genome}

        jasmine file_list={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.delly.simple_dupToIns_normalizeTypes.vcf max_dist=200 --allow_intrasample --comma_filelist --use_end --ignore_strand --nonlinear_dist out_file={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.delly.refined.vcf
        jasmine file_list={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.gridss.simple_dupToIns_normalizeTypes.vcf max_dist=200 --allow_intrasample --comma_filelist --use_end --ignore_strand --nonlinear_dist out_file={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.gridss.refined.vcf
        jasmine file_list={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.lumpy.simple_dupToIns_normalizeTypes.vcf max_dist=200 --allow_intrasample --comma_filelist --use_end --ignore_strand --nonlinear_dist out_file={params.out_dir}/{wildcards.tumor}.{wildcards.normal}.lumpy.refined.vcf

        ls {params.out_dir}/*.refined.vcf > {params.out_dir}/refined_files.txt
        jasmine file_list={params.out_dir}/refined_files.txt out_file={output.merged} --use_end --ignore_strand

        bcftools filter -i 'SUPP>1' {output.merged} > {output.overlap}
        """

def get_svs_merged_by_tool(wildcards):
    if wildcards.cell == "colo829":
        libs = ['B','C','D']
    else:
        libs = ['A','B','C']
    return [f"analysis/svs/jasmine_merge/{wildcards}/{lib}/{wildcards.cell}.{wildcards.cell}_bl.merged.supp2.vcf" for lib in libs]

rule jasmine_merge_libs:
    input:
        get_svs_merged_by_tool
    output:
        merged = "analysis/svs/jasmine_merge/{cell}/{cell}.final_merged.vcf",
        overlap = "analysis/svs/jasmine_merge/{cell}/{cell}.final_merged.supp2.vcf",
        dup_to_ins =  "analysis/svs/jasmine_merge/{cell}/{cell}.final_merged.supp2_dupToIns.vcf",
    params:
        out_dir = "analysis/svs/jasmine_merge/{cell}",
        jasmine_env = config['jasmine_env']
    shell:
        """
        module load conda-envs/base
        module load bcftools/1.19
        conda activate {params.jasmine_env}
        
        ls {input} > {params.out_dir}/final_merged.txt 
        jasmine file_list={params.out_dir}/final_merged.txt out_file={output.merged} --use_end --ignore_strand

        bcftools filter -i 'SUPP>1' {output.merged} > {output.overlap}

        jasmine --dup_to_ins --postprocess_only out_file={output.overlap} out_dir=analysis/svs/jasmine_merge/{wildcards.cell}
        """


