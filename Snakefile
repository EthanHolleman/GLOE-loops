import pandas as pd
from pathlib import Path

drip_samples_df = pd.read_table(
    'samples/DRIP_samples.csv', sep=',').set_index('sample_name', drop=False)

DRIP_SAMPLES = list(drip_samples_df['sample_name'])
GLOW_SAMPLES = [str(p.name) for p in Path('rawdata/GLOE-seq').iterdir()]
print(GLOW_SAMPLES)

rule all:
    input:
        expand('output/intersections/{sample}.bed', sample=DRIP_SAMPLES)

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
        GLOE_beds=expand('output/bed_sorted/GLOE/{gloe_bed}.sorted.bed', gloe_bed=GLOW_SAMPLES)
    output:
        'output/intersections/{sample}.bed'
    shell:'''
    bedtools intersect -wa -wb -a {input.DRIP_bed} -b {input.GLOE_beds} -sorted -filenames
    '''







# want to intersect each against all DRIP data??



