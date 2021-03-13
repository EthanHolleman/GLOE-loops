#!/bin/bash -l


mkdir -p Logs
conda activate snakemake
snakemake -j 12 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.cpus} --mem {cluster.mem} -J {cluster.name} -o {cluster.output} -e {cluster.output} --mail-type ALL --mail-user {cluster.email}" --latency-wait 120 -p -k $@
