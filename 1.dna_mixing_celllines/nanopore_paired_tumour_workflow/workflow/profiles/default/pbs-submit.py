#!/usr/bin/env python3
import sys
import subprocess

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]

job_properties = read_job_properties(jobscript)

ncpus=""
mem=""
walltime=""
ngpus=""

if 'threads' in job_properties:
      ncpus = f"ncpus={job_properties['threads']}"


if 'resources' in job_properties:
    resources = job_properties['resources']
    if 'gpu' in resources:
        ngpus = f"ngpus={resources['gpu']},"
    elif 'nvidia_gpu' in resources:
        ngpus = f"ngpus={resources['nvidia_gpu']},"

    if 'mem' in resources:
        mem = f"mem={resources['mem']}gb"
    
    if 'walltime' in resources:
        walltime = f"walltime={resources['walltime']}:00:00"
    
## whether use gpu nodes
queue = f"-q gpu" if len(ngpus) > 1 else ""


cmd = f"qsub {queue} -r n -l {ngpus}{ncpus},{mem},{walltime} {jobscript}"

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

res = res.stdout.decode()
print(res.strip())
