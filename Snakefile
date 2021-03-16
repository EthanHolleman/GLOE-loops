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
        expand('output/{sample}/indirect/{sample}.indirect.sorted.trim.bed', sample=GLOE_SAMPLES)


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
    conda:
        'envs/fastqc.yml'
    input:
        'rawdata/GLOE/{sample}.fastq'
    output:
        directory('output/fastqc/{sample}')
    shell:'''
    mkdir {output}
    fastqc {input} -o {output}
    '''


rule download_primer_file:
    output:
        'output/primers/TruSeq3-SE.fa'
    shell:'''
    curl https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa \
    -o {output}
    '''


rule trimmomatic:
    conda: 
        'envs/trimmomatic.yml'
    input:
        reads='rawdata/GLOE/{sample}.fastq.gz',
        primers='output/primers/TruSeq3-SE.fa'
    output:
        'output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz'
    threads: 16
    shell:'''
    trimmomatic SE -threads {threads} -phred33 {input.reads} {output} \
    ILLUMINACLIP:{input.primers}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
    '''


rule map_reads:  # are these reads paired end (mates)?
    conda: 
        'envs/bowtie2.yml'
    input:
        sample_reads='output/{sample}/trimmomatic/{sample}.trimmed.fastq.gz',
        bt_index='rawdata/bowtie2/hg19_index',
    output:
        'output/{sample}/mapped/{sample}.sam'
    threads: 12
    shell:'''
    bowtie2 -q --very-fast -x {input.bt_index}/hg19 -U {input.sample_reads} \
    -p {threads} -S {output}
    '''


rule sam_to_bam:
    conda: 
        'envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sam'
    output:
        'output/{sample}/mapped/{sample}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''


rule sort_bam:
    conda: 
        'envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.bam'
    output:
        'output/{sample}/mapped/{sample}.sorted.bam'
    params:
        sort='output/{sample}_temp'
    threads: 16
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule trim_bam_bad_alignments:
    conda: 
        'envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.bam'
    output:
        'output/{sample}/mapped/{sample}.trim.bam'
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    '''


rule sort_trimmed_bam:
    conda: 
        'envs/samtools.yml'
    input:
        temp('output/{sample}/mapped/{sample}.trim.bam')
    output:
        temp('output/{sample}/mapped/{sample}.sorted.trim.bam')
    threads: 16
    params:
        sort='output/alignment/{sample}_trim_temp'
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule flagstat_bam:
    conda: 
        'envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bam'
    output:
        'output/{sample}/mapped/{sample}.flagstat.sorted.trim.bam'
    shell:'''
    samtools flagstat {input} > {output}
    '''


rule index_sorted_bam:
    conda: 
        'envs/samtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bam'
    output:
        'output/{sample}/mapped/{sample}.sorted.trim.index'
    shell:'''
    samtools index {input} {output}
    '''

# Rules below are the bam2bedI module broken up

rule bam_to_bed:
    conda: 
        'envs/bedtools.yml'
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bam'
    output:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    shell:'''
    bedtools bamtobed -i {input} > {output}
    '''


rule indirect_mode:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    output:
        sites=temp('output/{sample}/indirect/{sample}.sites.indirect.trim.bed')
    shell:"""
    perl scripts/indirect_mode.pl {input} > {output.sites}
    """


rule sort_indirect_mode:
    input:
        'output/{sample}/indirect/{sample}.sites.indirect.trim.bed'
    output:
        temp('output/{sample}/indirect/{sample}.sites.indirect.sorted.trim.bed')
    shell:'''
    sort -k1,1 -k2,2n -k 6 {input} > {output}
    '''


# https://www.geeksforgeeks.org/awk-command-unixlinux-examples/
rule get_second_column_indirect:
    input:
        'output/{sample}/indirect/{sample}.sites.indirect.sorted.trim.bed'
    output:
        temp('output/{sample}/indirect/{sample}.sites.indirect.sorted.col2.trim.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule indirect_big_awk:
    input:
        'output/{sample}/indirect/{sample}.sites.indirect.sorted.col2.trim.bed'
    output:
        'output/{sample}/indirect/{sample}.indirect.sorted.trim.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """


rule seperate_forward_strand:
    input:
        'output/alignment/{sample}/indirect/{sample}.indirect.sorted.trim.bed'
    output:
        'output/alignment/{sample}/indirect/{sample}.fwd.indirect.sorted.trim.bed'
    
    shell:'''
    grep "+" {input} > {output}
    '''


rule seperate_reverse_strand:
    input:
        'output/alignment/{sample}/indirect/{sample}.indirect.sorted.trim.bed'
    output:
        'output/alignment/{sample}/indirect/{sample}.rev.indirect.sorted.trim.bed'
    shell:'''
    grep "-" {input} > {output}
    '''

# end bam2bedI rules 

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
    conda:
        'envs/bigWigToBedGraph.yml'
    input:
        'rawdata/DRIP/{sample}'
    output:
        'output/bed/DRIP/{sample}.bed'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule GLOE_bw_to_bed:
    conda:
        'envs/bigWigToBedGraph.yml'
    input:
        'rawdata/GLOE-seq/{sample}'
    output:
        'output/bed/GLOE/{sample}.bed'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule sort_DRIP_bed_files:
    conda:
        'envs/bedtools.yml'
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
    conda:
        'envs/bedtools.yml'
    input:
        DRIP_bed='output/trim_bed/DRIP/{sample}.trim.bed',
        GLOE_beds=expand('output/bed_sorted/GLOE/{gloe_bed}.sorted.bed', gloe_bed=GLOE_SAMPLES)
    output:
        'output/intersections/{sample}.bed'
    shell:'''
    bedtools intersect -wa -wb -a {input.DRIP_bed} -b {input.GLOE_beds} \
    -sorted -filenames > {output}
    '''




