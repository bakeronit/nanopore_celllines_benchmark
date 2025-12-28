for i in "" _10 _20 _30 _40 _50 _60 _70 _80 _90; do
    python gnomad_anno.py ../work/analysis/snvs/clairS/R10/sup/COLO829${i}.COLO829_BL/output.vcf.gz
    python gnomad_anno.py ../work/analysis/snvs/deepsomatic/R10/sup/COLO829${i}.COLO829_BL/output.somatic.vcf.gz
    python gnomad_anno.py ../work/analysis/snvs/clairS/R10/sup/HCC1937${i}.HCC1937_BL/output.vcf.gz
    python gnomad_anno.py ../work/analysis/snvs/deepsomatic/R10/sup/HCC1937${i}.HCC1937_BL/output.somatic.vcf.gz
done


for file in `ls ../../2.simulate_sequencing_depth/analysis/snvs/deepsomatic/*/output.somatic.vcf.gz`; do
for file in `ls ../../2.simulate_sequencing_depth/analysis/snvs/clairS/*/output.vcf.gz`; do
    python gnomad_anno.py $file
done

