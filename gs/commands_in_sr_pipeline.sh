#step 1 run GTAK
gatk --java-options "-Xmx8G" HaplotypeCaller \
  -R reference.fasta \
  -I sample.bam \
  -O sample.g.vcf \
  -L target_contig > gatk.log 2>&1

#step2. prepare and run qSNP
java -Xmx1G -cp qsnp.jar CreateQsnpIni \
  --ref reference.fasta \
  --mode somatic \
  --testBam tumor.bam \
  --testSample tumor \
  --controlBam normal.bam \
  --controlSample normal \
  --iniFile qsnp.ini

# Run qSNP
java -Xmx8G -jar qsnp.jar -i qsnp.ini -log qsnp.log

# Generate INI file
java -Xmx1G -cp q3indel.jar CreateQ3IndelIni \
  --ref reference.fasta \
  --mode somatic \
  --testBam tumor.bam \
  --testSample tumor \
  --controlBam normal.bam \
  --controlSample normal \
  --iniFile q3indel.ini

# Run Q3Indel
java -Xmx8G -jar q3indel.jar -i q3indel.ini -log q3indel.log

# Generate output filename
outputVcf="tumor_vs_normal.merged.vcf"

# Merge VCFs (simplified illustration)
cat snv.vcf indel.vcf | sort -k1,1 -k2,2n > $outputVcf

# Compress and index
bgzip $outputVcf
tabix -p vcf $outputVcf.gz

java -Xmx8G -jar qannotate.jar \
  --mode vcf2maf \
  -i tumor_vs_normal.merged.vcf.gz \
  -o filtered.vcf \
  --log filter.log \
  -d reference.fasta
  
## this will results in *sp.snpVcf.vcf.gz and *sp.indelVcf.vcf.gz as somatic passed variants.