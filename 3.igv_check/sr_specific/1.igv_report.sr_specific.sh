create_report colo829.sr_specific.bedpe \
    --genome hg38 --flanking 1000 \
    --tracks bam_files/LR/COLO829.bam bam_files/LR/COLO829_BL.bam bam_files/SR/COLO829.bam bam_files/SR/COLO829_BL.bam\
    --output colo829.sr_specific.igv.html


create_report hcc1937.sr_specific.bedpe \
    --genome hg38 --flanking 1000 \
    --tracks bam_files/LR/HCC1937.bam bam_files/LR/HCC1937_BL.bam bam_files/SR/HCC1937.bam bam_files/SR/HCC1937_BL.bam\
    --output hcc1937.sr_specific.igv.html