## already submitted jobs, change memory based on chromosome size
get_mem() {
    chrom=$1
    fai=~/working/data/genome/chm13/chm13v2.0.fa.fai
    base_chr="chr22"  # I did run with chr22 which used ~70GB memory
    base_gb=70
    size=$(awk -v c="$chrom" '$1==c{print $2}' "$fai")
    base_size=$(awk -v c="$base_chr" '$1==c{print $2}' "$fai")

    if [[ -z "$size" || -z "$base_size" ]]; then
        echo "chromosome not found in $fai" >&2
    fi
    mem=$(awk -v s="$size" -v sb="$base_size" -v gb="$base_gb" 'BEGIN {printf "%d\n", gb*s/sb}')
    min_mem=30
    if (( mem < min_mem )); then
        mem=$min_mem
    fi
    echo "$mem"
}

export -f get_mem
for jobid in $(qstat -u jiaZ|grep -w Q | awk '{print $1}'); do
    echo "processing jobid: $jobid"
    script_submitted=$(qstat -fx -F json $jobid | jq --arg id "$jobid" '.Jobs[$id]."Submit_arguments"' | awk '{print $NF}' | sed 's/"//')
    chrom=$(grep properties $script_submitted |  sed 's/^# properties = //' | jq '.wildcards.chrom'| sed 's/"//g')
    new_mem=$(get_mem $chrom)
    echo "change requested memory of $chrom running job $jobid to ${new_mem}G"
    #qalter -l mem=${new_mem}gb $jobid
    #echo "successfully changed memory to ${new_mem}G for job $jobid"
done
