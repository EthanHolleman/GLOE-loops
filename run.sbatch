#!/bin/bash -login
#SBATCH -p high                      # partition, or queue, to assign to
#SBATCH -J GLOE_SNAKE         # name for job
#SBATCH -o j%j.out # file to write STDOUT to
#SBATCH -e j%j.err # file to write error output to
#SBATCH -N 1                        # one "node", or computer
#SBATCH -n 1                        # one task for this node
#SBATCH -c 1                        # one cores per task
#SBATCH -t 48:00:00                  
#SBATCH --mem=10gb                  # ask for no more than 10 GB of memory
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etholleman@ucdavis.edu

conda activate snakemake
./run.sh  
# run the file that runs the snakemake workflow which
# submits more jobs to the cluster
