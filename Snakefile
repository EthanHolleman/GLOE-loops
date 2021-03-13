import pandas as pd
import os
from pathlib import Path

drip_samples_df = pd.read_table(
    'samples/DRIP_samples.csv', sep=',').set_index('sample_name', drop=False)

gloe_reads_df = pd.read_table(
    'samples/GLOE_samples.csv', sep=',').set_index('Sample Name', drop=False)

gloe_reads_df_url = gloe_reads_df[gloe_reads_df["url"].isnull() == False]
# only rows with SRA urls should be included

HG_LINK = 'https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/'
HG_NUMBERS = list(range(1, 23)) + ['X', 'Y']
DRIP_SAMPLES = list(drip_samples_df['sample_name'])
GLOE_SAMPLES = list(gloe_reads_df_url['Sample Name'])

hg_dict = {f'chr{chr_number}': os.path.join(HG_LINK, f'chr{chr_number}.fa.gz') 
           for chr_number in HG_NUMBERS}
hg_chromosomes = hg_dict.keys()


rule all:
    input:
        #expand('output/intersections/{sample}.bed', sample=DRIP_SAMPLES)
        #expand('rawdata/GLOE/{sample}.fastq', sample=GLOE_SAMPLES)
        expand('output/fastqc/{sample}', sample=GLOE_SAMPLES)


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
        downloaded_zip='rawdata/bowtie2/hg19_index/hg19.zip'
    shell:'''
    curl https://genome-idx.s3.amazonaws.com/bt/hg19.zip -o {output.downloaded_zip}
    unzip {output.downloaded_zip}
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

rule map_reads:  # are these reads paired end (mates)?
    input:
        sample_reads='rawdata/GLOE/{sample}.fastq',  # might not actually be this
        bw_index='rawdata/bowtie2/hg19_index',

    output:
        'output/aligment/{sample}.sam'
    shell:'''
    bowtie2 -x {input.bw_index} -s {output}
    '''

# rule expand_hg:
#     input:
#         expand('rawdata/hg/{chr_name}.fa.gz')


# rule download_hg:
#     output:
#         'rawdata/hg/{chr_name}.fa.gz'
#     params:
#         download_link = lambda wildcards: hg_dict[wildcards.chr_number]
#     shell:'''
#     curl -L {params.download_link} -o {output}
#     '''

# rule cat_hg_chromosomes:
#     input:
#         expand('rawdata/hg/{chr_name}.fa.gz', chr_name=hg_chromosomes)
#     output:
#         'rawdata/hg/complete.fa.gz'  # gz files can be directly concatenated
#     shell:'''
#     cat {input} > {output}
#     '''
    

rule expand_gloe_reads:
    input:
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
        'rawdata/GLOE/{sample}.fastq'
    shell:'''
    fastq-dump -Z {input} > {output}
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







# want to intersect each against all DRIP data??



