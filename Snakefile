import pandas as pd
from pathlib import Path

drip_samples_df = pd.read_table('samples/DRIP_samples.csv', sep=',').set_index('sample_name', drop=False)
DRIP_SAMPLES = list(drip_samples_df['sample_name'])
GLOW_SAMPLES = [str(p) for p in Path('rawdata/GLOE-seq').iterdir()]

rule all:
    input:
        expand('rawdata/DRIP/{sample}', sample=DRIP_SAMPLES)

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

rule bw_to_bed:
    input: '{sample}'
    output: 'output/bed/{sample}.bed'
    shell:'''
    bigWigtoBedGraph {input} {output}
    '''



# want to intersect each against all DRIP data??



