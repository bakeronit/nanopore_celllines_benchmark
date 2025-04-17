#!/usr/bin/env python3 -u
import subprocess
import sys
from cyvcf2 import VCF
from pathlib import Path


def get_snp_af_in_gnomad(mutation, queries):
    chrom, start, end, ref, alt = mutation.split('_')
    f_mutation = f"{chrom}_{start}_{ref}_{alt}" # they don't have a matched format: eg 1 -> chr1, chrom_start_end_ref_alt -> chrom_pos_ref_alt

    if queries[0] == '' or len([q for q in queries if q.split('\t')[0] == f_mutation]) == 0: # if queries[0] == '' is True the second statement won't process, so should have no error ab    out index out of range.
        af = 0
        return af
    query = [q for q in queries if q.split('\t')[0] == f_mutation][0]
    gnomad_mutation = query.split('\t')[0]
    af = float(query.split('\t')[-1]) if query.split('\t')[-1] != '.' else 0 # there are cases without AF_grpmax info, I checked the AF of it usually is 0 in different populations
    return af

def main():
    gnomad_dir="/reference/data/gnomAD/gnomad-public/4.1/vcf/genomes/"
    query_string = r"'%CHROM\_%POS\_%REF\_%ALT\t%ID\t%INFO/AF_grpmax'"
    vcf_file = Path(sys.argv[1])
    out_file = vcf_file.parent / "snv_gnomad_af_anno.tsv"
    print(out_file)
    out = open(out_file, 'wt')
    primary_chroms = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY","chrM"] 

    for variant in VCF(vcf_file):
        if not variant.CHROM in primary_chroms or variant.FILTER == "LowQual" or not variant.is_snp:
            continue
        end_pos = variant.POS + len(variant.REF)
        mutation = f"{variant.CHROM}_{variant.POS}_{end_pos}_{variant.REF}_{variant.ALT[0]}"
        region = f"{variant.CHROM}:{variant.POS}-{end_pos}"
        gnomad_vcf = gnomad_dir + f"gnomad.genomes.v4.1.sites.{variant.CHROM}.vcf.bgz"
        cmd = f"bcftools query -r {region} -f {query_string} {gnomad_vcf}"
        queries = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip('\n').split('\n')
        gnomad_af = get_snp_af_in_gnomad(mutation, queries)
        #af_col = "VAF" if "deepsomatic" in vcf_file.parts else "AF"
        #print(f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{variant.QUAL}\t{variant.format(af_col)[0][0]}\t{variant.format('DP')[0][0]}\t{gnomad_af}", file=out, flush=True)
        print(f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gnomad_af}", file=out, flush=True)
    out.close()

if __name__ == "__main__":
    main()
