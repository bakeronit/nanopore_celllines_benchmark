colo829_gs=/mnt/backedup/home/jiaZ/working/general/goldstandard/structural_variation/analysis/svs/jasmine_merge/colo829/colo829.final_merged.supp2.vcf
hcc1937_gs=/mnt/backedup/home/jiaZ/working/general/goldstandard/structural_variation/analysis/svs/jasmine_merge/hcc1937/hcc1937.final_merged.supp2.vcf

module load bcftools/1.19
bcftools view -i ID=@colo829.sr_specific.id $colo829_gs | \
    bcftools query -f "%CHROM\t%POS\t%POS\t%INFO/CHR2\t%INFO/END\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN" |\
    grep -v "#" > colo829.sr_specific.bedpe

bcftools view -i ID=@hcc1937.sr_specific.id $hcc1937_gs | \
    bcftools query -f "%CHROM\t%POS\t%POS\t%INFO/CHR2\t%INFO/END\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN" |\
    grep -v "#" > hcc1937.sr_specific.bedpe
