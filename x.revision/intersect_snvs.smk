"""
This analysis aim to evalute the recall and precision in ensembled calls as suggested by Reviewer #1
It is valuable to inspect how many false positive (precision) in the intersection of two tools, and how many false negative (recall) in the union calls of two tools.

The precision reflects how acurate we can be from a high-confidence calls.
The recall reflects how many calls will be missed from LRS regardless of tools.
"""

include: "helper.smk"

rule all:
    input:
        [f"analysis/snvs/intersections/{pair['tumour']}.{pair['normal']}/summary.txt" for pair in pairs]

rule bcftools_isec:
    input:
        clairs = purity_workdir / "analysis/snvs/clairS/R10/sup/{sample_t}.{sample_n}/output.vcf.gz",
        dpsomatic = purity_workdir / "analysis/snvs/deepsomatic/R10/sup/{sample_t}.{sample_n}/output.somatic.vcf.gz"
    output:
        "analysis/snvs/intersections/{sample_t}.{sample_n}/0000.vcf",
        "analysis/snvs/intersections/{sample_t}.{sample_n}/0001.vcf",
        "analysis/snvs/intersections/{sample_t}.{sample_n}/0002.vcf",
        "analysis/snvs/intersections/{sample_t}.{sample_n}/0003.vcf",
        "analysis/snvs/intersections/{sample_t}.{sample_n}/sites.txt",
        "analysis/snvs/intersections/{sample_t}.{sample_n}/README.txt",
    params:
        pdir = "analysis/snvs/intersections/{sample_t}.{sample_n}"
    envmodules:
        "bcftools/1.16"
    threads: 1
    resources:
        mem=2,
        walltime=1
    shell:
        """
        bcftools isec -p {params.pdir} {input.clairs} {input.dpsomatic} 
        """

rule intersect_summary:
    input:
        sites = "analysis/snvs/intersections/{sample_t}.{sample_n}/sites.txt",
        gs_sites = get_gs_sites
    output:
        "analysis/snvs/intersections/{sample_t}.{sample_n}/summary.txt"
    threads: 1
    resources:
        mem=1,
        walltime=1
    run:
        import gzip
        from pathlib import Path
        from collections import defaultdict
        gs_sites_file = Path(input.gs_sites)
        sites_file = Path(input.sites)
        sample_id = sites_file.parent.stem
        def valid_chrom(chrom: str) -> bool:
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            if chrom in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                return True
            return False

        def get_sites(sites_file: Path | str):
            _open = gzip.open if str(sites_file).endswith(".gz") else open
            sites = defaultdict(set)
            with _open(sites_file, 'rt') as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    chrom, pos, ref, alt, isec = line.strip().split("\t")
                    if len(ref) != 1 or len(alt) != 1 or not valid_chrom(chrom):
                        continue
                    site = f"{chrom}:{pos}:{ref}:{alt}"
                    sites[isec].add(site)
            return sites

        gs_sites = get_sites(gs_sites_file)["111"] | get_sites(gs_sites_file)["101"] | get_sites(gs_sites_file)["110"] | get_sites(gs_sites_file)["011"]
        sites = get_sites(sites_file)
        with open(output[0], 'wt') as out:
            out.write("\t".join(["sample", "intersection", "uniq_clairs", "uniq_dpsomatic", "union", "intersection_precision", "intersection_recall", "union_recall"]) + "\n")
            intersection = sites["11"]
            uniq_clairs = len(sites["10"])
            uniq_dpsomatic = len(sites["01"])
            union = intersection | sites["10"] | sites["01"]
            intersection_precision = len(intersection & gs_sites) / len(intersection) if len(intersection) > 0 else 0
            intersection_recall = len(intersection & gs_sites) / len(gs_sites)
            union_recall = len(union & gs_sites) / len(gs_sites)
            out.write("\t".join([sample_id, str(len(intersection)), str(uniq_clairs), str(uniq_dpsomatic), str(len(union)), f"{intersection_precision:.4f}", f"{intersection_recall:.4f}", f"{union_recall:.4f}"]) + "\n")
                
