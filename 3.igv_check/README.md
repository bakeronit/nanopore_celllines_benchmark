# SR-Specific Structural Variant IGV check

### Software Requirements
- `bcftools` v1.19
- [igv-reports](https://github.com/igvteam/igv-reports)

### Identification of SR-Specific SVs
We extracted SVs that were uniquely detected in SR data for both cell lines:
- COLO829: [List of SR-specific SVs](sr_specific/colo829.sr_specific.id)
- HCC1937: [List of SR-specific SVs](sr_specific/hcc1937.sr_specific.id)

#### 1. Extraction of SR-Specific SVs as BEDPE Files

The following commands were used to extract SR-specific SVs from the gold standard VCF files and convert them to BEDPE format:

```bash
# Extract COLO829 SR-specific SVs as BEDPE
bcftools view -i ID=@colo829.sr_specific.id $colo829_gs | \
    bcftools query -f "%CHROM\t%POS\t%POS\t%INFO/CHR2\t%INFO/END\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN" | \
    grep -v "#" > colo829.sr_specific.bedpe

# Extract HCC1937 SR-specific SVs as BEDPE
bcftools view -i ID=@hcc1937.sr_specific.id $hcc1937_gs | \
    bcftools query -f "%CHROM\t%POS\t%POS\t%INFO/CHR2\t%INFO/END\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN" | \
    grep -v "#" > hcc1937.sr_specific.bedpe
```

#### 2. Visualization with IGV-Reports

To visually inspect these SR-specific SVs, we generated interactive HTML reports using `igv-reports`:

```bash
# Generate IGV report for COLO829 SR-specific SVs
create_report colo829.sr_specific.bedpe \
    --genome hg38 --flanking 1000 \
    --tracks bam_files/LR/COLO829.bam bam_files/LR/COLO829_BL.bam bam_files/SR/COLO829.bam bam_files/SR/COLO829_BL.bam \
    --output colo829.sr_specific.igv.html

# Generate IGV report for HCC1937 SR-specific SVs
create_report hcc1937.sr_specific.bedpe \
    --genome hg38 --flanking 1000 \
    --tracks bam_files/LR/HCC1937.bam bam_files/LR/HCC1937_BL.bam bam_files/SR/HCC1937.bam bam_files/SR/HCC1937_BL.bam \
    --output hcc1937.sr_specific.igv.html
```

## Results

### Visual Inspection
The generated IGV reports allow for a visual examination of each SR-specific SV, and we found them might be false positive in SR.
e.g. [COLO829 IGV Report](sr_specific/colo829.sr_specific.igv.html). 
*HTML file size for HCC1937 are too big to push to github.*

