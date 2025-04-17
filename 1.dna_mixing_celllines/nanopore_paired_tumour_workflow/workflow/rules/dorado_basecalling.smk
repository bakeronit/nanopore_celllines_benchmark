DORADO = config['dorado_latest']

rule dorado_basecalling:
    input:
        unpack(get_basecalling_models),
        pod5 = "raw_pod5/{flowcell}/{sample}/{run}"
    output:
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam"
    benchmark:
        "benchmarks/dorado/{flowcell}.{sample}.{run}.{mode}.benchmark.txt"
    log:
        "logs/dorado/{flowcell}.{sample}.{run}.{mode}.log"
    threads: 8
    resources:
        nvidia_gpu = 1,
        walltime = 200,
        mem = 16
    envmodules:
        "samtools/1.17"
    shell:
        """
        {DORADO} basecaller {input.cn_model} {input.pod5} \
        --modified-bases-models {input.remora_model} \
        --recursive --device "cuda:all"  2>{log} | samtools view -b -o {output} -
        """