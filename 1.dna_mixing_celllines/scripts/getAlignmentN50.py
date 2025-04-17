import pysam
import sys
from numpy import cumsum, argmax
from pathlib import Path

def getN50(num_list):
    sorted_list = sorted(num_list, reverse=True)
    csum = cumsum(sorted_list)
    return sorted_list[argmax(csum >= csum[-1] * 0.5)] 

bam_file_path = Path(sys.argv[1])

alignment_length = []
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    for read in bam_file:
        alignment_length.append(read.query_alignment_length)

N50 = getN50(alignment_length)

sample = bam_file_path.name.removesuffix(".bam")
print(f"{sample}\t{N50}")
