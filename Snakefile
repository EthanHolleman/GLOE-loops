import pandas as pd
from pathlib import Path

drip_samples_df = pd.read_table(
    'samples/DRIP_samples.csv', sep=',').set_index('sample_name', drop=False)
DRIP_SAMPLES = list(drip_samples_df['sample_name'])
GLOW_SAMPLES = [str(p) for p in Path('rawdata/GLOE-seq').iterdir()]

rule all:
    input:
        expand('output/trim_bed/DRIP/{sample}', sample=DRIP_SAMPLES)

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

use rule DRIP_bw_to_bed as GLOE_bw_to_bed with:
    input:
        'rawdata/GLOE-seq/{sample}'
    output:
        'output/bed/GLOE/{sample}.bed'

# want to intersect all DRIP against all GLOE

rule trim_DRIP_peaks:
    input:
        'output/bed/DRIP/{sample}.bed'
    output:
        'output/trim_bed/DRIP/{sample}'
    shell:'''
    Rscript scripts/peak_sep.R {input} {output}
    '''

# rule test_input:
#     input:
#         DRIP=expand('output/bed/DRIP/{drip_sample}.bed', sample=DRIP_SAMPLES),
#         GLOE=expand('output/bed/GLOE/{gloe_sample}.bed', sample=DRIP_SAMPLES)

# rule intersect_all:
#     output:
#         'output/intersections/{sample_a}-{sample_b}'
#     shell:'''
#     bedtools intersect -wa -wb -a {drip_sample} -b {GLOE} -sorted -filenames
#     '''







# want to intersect each against all DRIP data??



