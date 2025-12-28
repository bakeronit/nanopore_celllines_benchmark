from pathlib import Path
import pandas as pd

purity_workdir = Path("../1.dna_mixing_celllines/work")
depth_workdir =  Path("../2.simulate_sequencing_depth")

config_dir = Path("../1.dna_mixing_celllines/nanopore_paired_tumour_workflow/config") 
data_list = [pd.read_csv(f) for f in config_dir.glob("sample*.csv")]
samples_df = pd.concat(data_list, ignore_index=True)

wildcard_constraints:
    sample = "|".join(samples_df['sample_id'].unique()),
    sample_t = "|".join(samples_df[samples_df['type'] == 'tumour']['sample_id'].unique()),
    sample_n = "|".join(samples_df[samples_df['type'] == 'normal']['sample_id'].unique()),
    depth = "|".join(["60x","45x","30x","15x"]),
    variant_type = "|".join(["snvs","indels"])

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

def get_gs_sites(wildcards):
    cellline = wildcards.sample_t.split("_")[0].lower()
    return f"../gs/vcfs/{cellline}/isec_hom100_snv_dir/sites.txt.gz"

pairs = generate_paired_samples(samples_df)

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
