import pandas as pd
from pathlib import Path
import sys

samples_df = pd.read_csv(config['samples'])

wildcard_constraints:
    sample="|".join(samples_df["sample_id"].unique()),
    run="|".join(samples_df["flowcell_id"].unique())

MODE = config['basecalling_mode']

def get_basecalling_output():
    output = expand(("analysis/ubam/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '/' + samples_df['flowcell_id'] + '.ubam'), mode=MODE)
    return output

def get_alignment_output():
    return expand(("analysis/bam/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '.bam').unique(), mode=MODE)

def generate_paired_samples(df):
    """generate tumour-normal pairs for each donor"""
    donor_flowcell_df = df[['donor_id','flowcell_version']].drop_duplicates() # for each donor, and potentially each flowcell version
    pairs = []
    for index, row in donor_flowcell_df.iterrows():
        donor_id = row['donor_id']
        flowcell = row['flowcell_version']
        tumour_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'tumour')]['sample_id'].unique().tolist()
        normal_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'normal')]['sample_id'].unique().tolist()
        if len(tumour_sample) > 1 and len(normal_sample) > 1:
            raise ValueError('Multiple tumour vs multiple normal samples for one donor, could not determine the paired samples')
            sys.exit()
        pairs += [{'donor': donor_id, 'tumour': t, 'normal': n, 'flowcell_version': flowcell} for t in tumour_sample for n in normal_sample]
    return pairs

def get_snv_calling_output(tool):
    if tool in ['clairS','deepsomatic']:
        pairs = generate_paired_samples(samples_df)
        output = [f"analysis/snvs/{tool}/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in pairs for m in MODE]
    if tool == 'clair3':
        clair3_path = Path('analysis/snvs/clair3')
        normal_samples_df = samples_df[samples_df['type']=='normal']
        output = expand(( clair3_path/normal_samples_df['flowcell_version']/'{mode}'/normal_samples_df['sample_id']/'phased_merge_output.vcf.gz').unique(), mode=MODE)
    if tool == 'deepvariant':
        deepvariant_path = Path('analysis/snvs/deepvariant/R10')  # deepvariant only works with R10 data.
        normal_samples = samples_df[(samples_df['flowcell_version']=='R10') & (samples_df['type']=='normal')]['sample_id']
        output = expand(( deepvariant_path/'{mode}'/normal_samples/(normal_samples + ".vcf.gz")).unique(), mode=MODE)
        output += expand(( deepvariant_path/'{mode}'/normal_samples/(normal_samples + ".g.vcf.gz")).unique(), mode=MODE)
    if tool == 'pepper':
        pepper_path = Path('analysis/snvs/pepper')
        normal_samples_df = samples_df[samples_df['type']=='normal']
        output = expand(( pepper_path/normal_samples_df['flowcell_version']/'{mode}'/normal_samples_df['sample_id']/(normal_samples_df['sample_id'] + '.phased.vcf.gz')).unique(), mode=MODE) 
    return output

def get_sv_calling_output(tool):
    pairs = generate_paired_samples(samples_df)
    match tool:
        case "severus":
            output = [f"analysis/svs/severus/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/somatic_SVs/severus_somatic.vcf" for pair in pairs for m in MODE]
        case "savana":
            output = [f"analysis/svs/savana/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.classified.somatic.vcf" for pair in pairs for m in MODE]
        case "savana1.20":
            output = [f"analysis/svs/savana1.20/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.classified.somatic.vcf" for pair in pairs for m in MODE]
        case "nanomonsv":
            output = [f"analysis/svs/nanomonsv/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.nanomonsv.result.filt.pass.svtype.txt" for pair in pairs for m in MODE] 
            #output += [f"analysis/svs/nanomonsv/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.nanomonsv.sbnd.annot.proc.result.pass.txt" for pair in pairs for m in MODE]
        case "delly":
            output =  [f"analysis/svs/delly/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.vcf" for pair in pairs for m in MODE]
        case _:
            raise ValueError("Tool can only be severus, savana, delly and nanomonsv, check the config file")
            return None
    return output

def get_final_output():
    ## qc
    #qc_output = expand(("analysis/qc/basecalling/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '/' + samples_df['flowcell_id'] + '.summary.txt'), mode=MODE)
    qc_output = expand(("analysis/qc/bam/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '{suffix}'), mode=MODE, suffix=['.mosdepth.summary.txt','.mosdepth.global.dist.txt','.bamcov.txt','.bamN50.txt']) 

    ## snv calling
    snv_output = []
    for tool in config['germline_snv_calling']:
        snv_output += get_snv_calling_output(tool)
    for tool in config['somatic_snv_calling']:
        snv_output += get_snv_calling_output(tool)
    
    ## methylation
    mod_output = expand( ('analysis/mod/' + samples_df['flowcell_version'] + '/{mode}/' + samples_df['sample_id'] + '.bed.gz').unique(), mode=MODE)
    mod_output += expand( ('analysis/mod/' + samples_df['flowcell_version'] + '/{mode}/' + samples_df['sample_id'] + '.mod_summary.txt').unique(), mode=MODE)

    ## sv output
    pairs = generate_paired_samples(samples_df)
    sv_output = []
    for tool in ['severus','savana','nanomonsv','delly']:
        sv_output += get_sv_calling_output(tool)
    
    ## cnv output
    cnv_output = [f"analysis/cnvs/savana1.20/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}_fitted_purity_ploidy.tsv" for pair in pairs for m in MODE]
    final_output = mod_output + snv_output + sv_output + qc_output + cnv_output

    step = config['step']
    match step:
        case 'basecalling':
            return get_basecalling_output()
        case 'alignment':
            return get_alignment_output() + qc_output
        case 'snv':
            return snv_output + mod_output + qc_output
        case 'all':
            return final_output
        case _:
            raise ValueError("step can only be basecalling, alignment, snv, and all.")
            return None

            
## get both cannonical and mod model files for basecalling
def get_basecalling_models(wildcards):
    if wildcards.flowcell == "R9":
        if wildcards.mode == "hac":
            cn_model = config['model']['r9']['hac']['canonical']
            remora_model = config['model']['r9']['hac']['remora']
        elif wildcards.mode == "sup":
            cn_model = config['model']['r9']['sup']['canonical']
            remora_model = config['model']['r9']['sup']['remora']
    elif wildcards.flowcell == "R10":
        if wildcards.mode == "hac":
            cn_model = config['model']['r10']['hac']['canonical']
            remora_model = config['model']['r10']['hac']['remora']
        elif wildcards.mode == "sup":
            cn_model = config['model']['r10']['sup']['canonical']
            remora_model = config['model']['r10']['sup']['remora']
    else:
        print("Error: flowcell version has to be R9 or R10. mode has to be hac or sup")
        exit(0)
    
    return {"cn_model": cn_model, "remora_model": remora_model}

## get the path to raw output from sequencer, either fast5 or pod5
def get_raw_path(wildcards):

    raw_path = list(samples_df[ (samples_df['flowcell_id'] == wildcards.run) ].raw_path)
    
    return raw_path
    
## get all bam files of a sample
def get_bam_of_runs(wildcards):
    sample = wildcards.sample
    flowcell = wildcards.flowcell
    mode = wildcards.mode

    align_results_path = Path(f'analysis/bam/{flowcell}/{mode}/{sample}')
    flowcell_ids = set(samples_df[ (samples_df['flowcell_version'] == flowcell) & (samples_df['sample_id'] == sample) ].flowcell_id)

    return [align_results_path / f'{id}.bam' for id in flowcell_ids]

## get model file for clair3 snp calling
def get_clair3_model(wildcards):
    if wildcards.flowcell == 'R9':
        clair3_model = config['clair3']['model_r9']
    elif wildcards.flowcell == 'R10':
        clair3_model = config['clair3']['model_r10']
    else:
        print("Error: flowcell version has to be R9 or R10.")
        exit(0)
    
    return clair3_model

def get_phased_vcf(sample, flowcell, mode):
    donor_id = samples_df[samples_df['sample_id']==sample].donor_id.tolist()[0]
    normal_sample_id = samples_df[(samples_df['donor_id']==donor_id) & (samples_df['type']=='normal')].sample_id.tolist()[0]
    match config['phased_snv_from']:
        case 'clair3':
            return f"analysis/snvs/clair3/{flowcell}/{mode}/{normal_sample_id}/phased_merge_output.vcf.gz"
        case 'pepper':
            return f"analysis/snvs/pepper/{flowcell}/{mode}/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case 'deepvariant':
            return f"analysis/snvs/deepvariant/{flowcell}/{mode}/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case _:
            raise ValueError("haplotagged bam should only be generated from clair3, pepper, or deepvariant!")
            return None
