## this just create soft links for all pod5 files for each sample instead of copying them.
##!!!! if starts from fast5, you need to request more threads and resources.
## currently I change to 1 thread and less mem, walltime to avoid over requesting.
########
rule gather_rawdata:
    input:
        get_raw_path
    output:
        pod5_dir = directory("raw_pod5/{flowcell}/{sample}/{run}"),
        fin = "raw_pod5/{flowcell}/{sample}/{run}.done"
    threads: 1
    resources:
        walltime = 1,
        mem = 1
    shell:
        """
        if [[  -z `find {input} -type f -name *.pod5` ]]; then
            module load longread/python39venv
            mkdir -p {output.pod5_dir}
            pod5 convert fast5 --recursive --strict --threads {threads} --output {output.pod5_dir} {input}
            touch {output.fin}
        else
            mkdir -p {output.pod5_dir}
            find {input} -type f -name *.pod5 | xargs -I {{}} ln -s {{}} {output.pod5_dir}
            touch {output.fin}
        fi
        """
