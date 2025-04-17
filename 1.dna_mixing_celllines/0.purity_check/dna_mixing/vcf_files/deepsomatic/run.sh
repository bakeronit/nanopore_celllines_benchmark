module load bcftools/1.19 htslib/1.19.1

for i in `ls /mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/*/output.somatic.vcf.gz`;do 
    sample=`echo $i|cut -d"/" -f16|cut -d"." -f1`
    #ln -s $i ${sample}.vcf.gz; ln -s ${i}.tbi ${sample}.vcf.gz.tbi
    bcftools view --types snps $i |bgzip > ${sample}.vcf.gz
    tabix -p vcf ${sample}.vcf.gz
done

mkdir -p cn{1,2}

for i in COLO829*.vcf.gz;do 
    sample=`basename $i .vcf.gz`
    bcftools view -R ../../copynumber/colo829.cn2.tsv $i > cn2/$sample.vcf
    bcftools view -R ../../copynumber/colo829.cn1.tsv $i > cn1/$sample.vcf
done 

for i in HCC1937*.vcf.gz;do 
    sample=`basename $i .vcf.gz`
    bcftools view -R ../../copynumber/hcc1937.cn2.tsv $i > cn2/$sample.vcf
    bcftools view -R ../../copynumber/hcc1937.cn1.tsv $i > cn1/$sample.vcf
done 
