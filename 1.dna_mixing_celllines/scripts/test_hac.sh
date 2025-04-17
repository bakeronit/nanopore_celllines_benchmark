#PBS -N bc_hac
#PBS -r n
#PBS -q gpu
#PBS -l ngpus=1,ncpus=8,mem=30gb,walltime=40:00:00
#PBS -m ae

cd /mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/scripts
 
dorado=/mnt/backedup/home/jiaZ/working/local/dorado/dorado-0.5.1-linux-x64/bin/dorado 
workdir=/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work
module load samtools/1.17

#mkdir -p ${workdir}/analysis/ubam/R10/hac/COLO829
#$dorado basecaller /mnt/backedup/home/jiaZ/working/data/ont_models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
#    ${workdir}/raw_pod5/R10/COLO829/PAQ39011 \
#    --modified-bases-models /mnt/backedup/home/jiaZ/working/data/ont_models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0_5mCG_5hmCG@v1 \
#     --recursive --device "cuda:all"  2> ${workdir}/logs/dorado/R10.COLO829.PAQ39011.hac.log | \
#     samtools view -b -o ${workdir}/analysis/ubam/R10/hac/COLO829/PAQ39011.ubam -

time $dorado basecaller /mnt/backedup/home/jiaZ/working/data/ont_models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
    ${workdir}/raw_pod5/R10/COLO829/PAQ22128 \
    --modified-bases-models /mnt/backedup/home/jiaZ/working/data/ont_models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0_5mCG_5hmCG@v1 \
     --recursive --device "cuda:all"  2> ${workdir}/logs/dorado/R10.COLO829.PAQ22128.hac.log | \
     samtools view -b -o ${workdir}/analysis/ubam/R10/hac/COLO829/PAQ22128.ubam -
