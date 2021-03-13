#!/bin/bash -l

module load bio3

mkdir -p Logs
snakemake -j 4 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.cpus} -J {cluster.name} -o {cluster.output} -e {cluster.output}"

