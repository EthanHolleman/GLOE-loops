import pandas as pd
import os
from pathlib import Path

drip_samples_df = pd.read_table(
    'samples/DRIP_samples.csv', sep=',').set_index('sample_name', drop=False)

gloe_reads_df = pd.read_table(
    'samples/GLOE_samples.csv', sep=',').set_index('Sample Name', drop=False)

gloe_reads_df_url = gloe_reads_df[gloe_reads_df["url"].isnull() == False]
# only rows with SRA urls should be included

DRIP_SAMPLES = list(drip_samples_df['sample_name'])
GLOE_SAMPLES = list(gloe_reads_df_url['Sample Name'])


rule all:
    input:
        #expand('output/intersections/{sample}.bed', sample=DRIP_SAMPLES)
        #expand('rawdata/GLOE/{sample}.fastq', sample=GLOE_SAMPLES)
        #expand('output/fastqc/{sample}', sample=GLOE_SAMPLES)
        #expand('output/alignment/{sample}.sam', sample=GLOE_SAMPLES)
        expand('output/alignment/{sample}/{sample}.sorted.trim.bam', sample=GLOE_SAMPLES)


rule expand_drip:
    input:
        expand('rawdata/DRIP/{sample}', sample=DRIP_SAMPLES)


rule download_drip:
    output:
        'rawdata/DRIP/{sample}'
    params:
        download_link = lambda wildcards: drip_samples_df.loc[wildcards.sample]['url']
    shell:'''
    curl -L {params.download_link} -o {output}
    '''


rule download_hg_bt_index:
    output:
        index_dir=directory('rawdata/bowtie2/hg19_index'),
        downloaded_zip='rawdata/bowtie2/hg19_index/hg19.zip',
    shell:'''
    mkdir --parents {output.index_dir}
    curl https://genome-idx.s3.amazonaws.com/bt/hg19.zip -o {output.downloaded_zip}
    unzip {output.downloaded_zip} -d {output.index_dir}
    '''


rule fastqc:
    input:
        'rawdata/GLOE/{sample}.fastq'
    output:
        directory('output/fastqc/{sample}')
    shell:'''
    mkdir {output}
    fastqc {input} -o {output}
    '''

rule trimmomatic:
    input:
        'rawdata/GLOE/{sample}.fastq.gz'
    output:
        'output/{sample}/{sample}.trimmed.fastq.gz'
    shell:'''
    trimmomatic SE -phred33 {input} {output} \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10, SLIDINGWINDOW:4:15, MINLEN:36
    '''


rule map_reads:  # are these reads paired end (mates)?
    input:
        sample_reads='output/{sample}/{sample}.trimmed.fastq.gz',
        bt_index='rawdata/bowtie2/hg19_index',
    output:
        'output/alignment/{sample}/{sample}.sam'
    threads: 12
    shell:'''
    bowtie2 -q --very-fast -x {input.bt_index}/hg19 -U {input.sample_reads} \
    -p {threads} -S {output}
    '''


rule sam_to_bam:
    input:
        'output/alignment/{sample}/{sample}.sam'
    output:
        'output/alignment/{sample}/{sample}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''


rule sort_bam:
    input:
        'output/alignment/{sample}/{sample}.bam'
    output:
        'output/alignment/{sample}/{sample}.sorted.bam'
    params:
        temp='output/alignment/{sample}_temp'
    threads: 16
    shell:'''
    mkdir --parents {params.temp}
    samtools sort -O bam -T {params.temp} --threads {threads} -o {output} {input}
    '''


rule trim_bam_bad_alignments:
    input:
        'output/alignment/{sample}/{sample}.sorted.bam'
    output:
        'output/alignment/{sample}/{sample}.trim.bam'
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    rm {input}
    '''


rule sort_trimmed_bam:
    input:
        'output/alignment/{sample}/{sample}.trim.bam'
    output:
        'output/alignment/{sample}/{sample}.sorted.trim.bam'
    threads: 16
    params:
        temp='output/alignment/{sample}_trim_temp'
    shell:'''
    mkdir --parents {params.temp}
    samtools sort -O bam -T {params.temp} --threads {threads} -o {output} {input}
    rm {input}
    '''


rule expand_gloe_reads:
    input:
        # somehow wildcard sample = sample_name.trimmed ??
        # where is trimmed coming from?
        expand('rawdata/GLOE/{sample}.sra', sample=GLOE_SAMPLES)


rule download_gloe_reads:
    output:
        'rawdata/GLOE/{sample}.sra'
    params:
        download_link = lambda wildcards: gloe_reads_df_url.loc[wildcards.sample]['url']
    shell:'''
    curl -L {params.download_link} -o {output}
    '''


rule dump_gloe_fastq:
    input:
        'rawdata/GLOE/{sample}.sra'
    output:
        'rawdata/GLOE/{sample}.fastq.gz'
    shell:'''
    fastq-dump -Z {input} | gzip > {output}
    '''


rule DRIP_bw_to_bed:
    input:
        'rawdata/DRIP/{sample}'
    output:
        'output/bed/DRIP/{sample}.bed'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule GLOE_bw_to_bed:
    input:
        'rawdata/GLOE-seq/{sample}'
    output:
        'output/bed/GLOE/{sample}.bed'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule sort_DRIP_bed_files:
    input:
        'output/bed/DRIP/{sample}.bed'
    output:
        'output/bed_sorted/DRIP/{sample}.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


use rule sort_DRIP_bed_files as sort_GLOE_bed_files with:
    input:
        'output/bed/GLOE/{sample}.bed'
    output:
         'output/bed_sorted/GLOE/{sample}.sorted.bed'
    

rule trim_DRIP_peaks:
    input:
        'output/bed_sorted/DRIP/{sample}.sorted.bed'
    output:
        'output/trim_bed/DRIP/{sample}.trim.bed'
    shell:'''
    Rscript scripts/peak_sep.R {input} {output}
    '''


rule intersect_all:
    input:
        DRIP_bed='output/trim_bed/DRIP/{sample}.trim.bed',
        GLOE_beds=expand('output/bed_sorted/GLOE/{gloe_bed}.sorted.bed', gloe_bed=GLOE_SAMPLES)
    output:
        'output/intersections/{sample}.bed'
    shell:'''
    bedtools intersect -wa -wb -a {input.DRIP_bed} -b {input.GLOE_beds} \
    -sorted -filenames > {output}
    '''




