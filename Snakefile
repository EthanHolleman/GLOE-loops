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
FWD_DRIP_SAMPLES = [sample for sample in DRIP_SAMPLES if '_pos' in sample]
REV_DRIP_SAMPLES = [sample for sample in DRIP_SAMPLES if '_neg' in sample]
GLOE_SAMPLES = list(gloe_reads_df_url['Sample Name'])

STRANDS = ['fwd', 'rev']
MODES = ['direct', 'indirect']
CLOSEST_MODES = ['all', 'first']



rule all:
    input:
        expand(
            'output/{sample}/random/{mode}.{sample}.{strand}.random.closest.first.bed', 
            sample=GLOE_SAMPLES, mode=MODES, strand=STRANDS
        ),
        expand(
            'output/plots/footloop_dist_vs_random/{mode}.{sample}.footloop_{footloop_strand}.gloe_{gloe_strand}.nick_dist.png',
            sample=GLOE_SAMPLES, mode=MODES, footloop_strand=STRANDS, gloe_strand=STRANDS
        ),
        expand(
            'output/footloop_cons/footloop_all.{strand}.con.bed',
            strand=STRANDS
        ),
        expand(
            'output/{sample}/footloop_con/{mode}.footloop_all.{strand}.{sample}.{gloe_strand}.con.sorted.first.closest.png',
            sample=GLOE_SAMPLES, mode=MODES, gloe_strand=STRANDS, strand=STRANDS
        ),

         

rule download_footloop_all:
    output:
        'rawdata/footloop/footloop_all.bed'
    shell:'''
    curl -L "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1061960477_AzLo24AGPf9cdZa26GHovMN4XI9v&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_ct_footLoopPeakALL_41&hgta_ctDesc=table+browser+query+on+ct_footLoopPeakALL_41&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED" -o {output}
    '''


rule seperate_footloop_strands:
    input:
        'rawdata/footloop/footloop_all.bed'
    output:
        fwd='rawdata/footloop/footloop_all.fwd.bed',
        rev='rawdata/footloop/footloop_all.rev.bed'
    shell:'''
    grep "+" {input} > {output.fwd}
    grep "-" {input} > {output.rev}
    '''


rule sort_footloop_strands:
    input:
        'rawdata/footloop/footloop_all.{direction}.bed'
    output:
        'rawdata/footloop/footloop_all.{direction}.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule window_footloop:
    input:
        footloop='rawdata/footloop/footloop_all.{direction_1}.sorted.bed',
        gloe='output/{sample}/{mode}/{sample}.{direction_2}.{mode}.sorted.trim.bed'
    output:
        'output/footloop_window/{sample}.{direction_2}.{direction_1}.window.{mode}.bed'
    params:
        w="1000"  # window size up and downstream
    shell:'''
    mkdir -p output/footloop_window/
    bedtools window -w {params.w} -a {input.footloop} -b {input.gloe} > {output}
    '''


rule plot_window_footloop:
    input:
        'output/footloop_window/{sample}.{direction_2}.{direction_1}.window.{mode}.bed'
    output:
        'output/plots/footloop_window/{sample}.{direction_2}.{direction_1}.png'
    shell:'''
    mkdir -p output/plots/footloop_window/
    /home/ethollem/anaconda3/bin/Rscript scripts/footloop_nick_dist.R {input} {output}
    '''


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


rule truncate_footloops:
    input:
        'rawdata/footloop/footloop_all.{strand}.bed'
    output:
        temp('output/truncated_footloop/footloop_all.{strand}.trunc.bed')
    shell:'''
    mkdir -p output/truncated_footloop/
    /home/ethollem/anaconda3/bin/Rscript scripts/truncate_footloops.R {input} {output}
    '''


rule sort_truncated_footloops:
    input:
        'output/truncated_footloop/footloop_all.{strand}.trunc.bed'
    output:
        'output/truncated_footloop/footloop_all.{strand}.trunc.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule closest_truncated_breaks_first:
    input:
        footloop='output/truncated_footloop/footloop_all.{footloop_strand}.trunc.sorted.bed',
        gloe= 'output/{sample}/{mode}/{sample}.{gloe_strand}.{mode}.sorted.trim.bed'
    output:
        'output/truncated_closests/footloop_all.{mode}.footloop_{footloop_strand}.gloe_{gloe_strand}.trunc.closest.first.{sample}.bed'
    shell:'''
    mkdir -p output/truncated_closests
    bedtools closest -t first -D ref -a {input.footloop} -b {input.gloe} > {output}
    '''

rule closest_truncated_breaks_all:
    # closest to the R-loop initiation site, in case of ties report them all
    # this is basically allow for counting the number of reads at the closest
    # location as ties would be of the same location
    input:
        footloop='output/truncated_footloop/footloop_all.{footloop_strand}.trunc.sorted.bed',
        gloe= 'output/{sample}/{mode}/{sample}.{gloe_strand}.{mode}.sorted.trim.bed'
    output:
        'output/truncated_closests/footloop_all.{mode}.footloop_{footloop_strand}.gloe_{gloe_strand}.trunc.closest.all.{sample}.bed'
    shell:'''
    mkdir -p output/truncated_closests
    bedtools closest -D ref -a {input.footloop} -b {input.gloe} > {output}
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
    bowtie2 -q -x {input.bt_index}/hg19 -U {input.sample_reads} \
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
        'output/{sample}/mapped/{sample}.trim.bam'
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


rule GLOW_bed_to_bam:
    input:
        bed='output/{sample}/{mode}/{sample}.{strand}.{mode}.sorted.trim.bed',
        genome='rawdata/hg19/hg19.chrom.sizes'
    output:
        temp('output/{sample}/{mode}/{sample}.{strand}.{mode}.trim.bam')
    shell:'''
    bedtools bedtobam -i {input.bed} -g {input.genome} > {output}
    '''

rule sort_GLOW_bam:
    input:
        'output/{sample}/{mode}/{sample}.{strand}.{mode}.trim.bam'
    output:
        'output/{sample}/{mode}/{sample}.{strand}.{mode}.trim.sorted.bam'
    shell:'''
    samtools sort {input} > {output}
    '''


rule GLOW_read_depth:
    conda:
        'envs/samtools.yml'
    input:
        'output/{sample}/{mode}/{sample}.{strand}.{mode}.trim.sorted.bam'
    output:
        'output/{sample}/coverage/{sample}.{strand}.{mode}.coverage.sorted.trim.bed'
    params:
        output_dir='output/{sample}/coverage'
    shell:'''
    mkdir -p {params.output_dir}
    samtools depth {input} > {output}
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

rule direct_mode:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    output:
        sites='output/{sample}/direct/{sample}.sites.direct.trim.bed'
    params:
        index_dir='output/{sample}/direct'
    shell:'''
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/direct_mode.pl {input} > {output.sites}
    '''

rule indirect_mode:
    input:
        'output/{sample}/mapped/{sample}.sorted.trim.bed'
    output:
        sites='output/{sample}/indirect/{sample}.sites.indirect.trim.bed'
    params:
        index_dir='output/{sample}/indirect'
    shell:"""
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/indirect_mode.pl {input} > {output.sites}
    """


rule sort_perl_mode_output:
    input:
        'output/{sample}/{mode}/{sample}.sites.{mode}.trim.bed'
    output:
        temp('output/{sample}/{mode}/{sample}.sites.{mode}.sorted.trim.bed')
    shell:'''
    sort -k1,1 -k2,2n -k 6 {input} > {output}
    '''


# https://www.geeksforgeeks.org/awk-command-unixlinux-examples/
rule get_second_column:
    input:
        'output/{sample}/{mode}/{sample}.sites.{mode}.sorted.trim.bed'
    output:
        temp('output/{sample}/{mode}/{sample}.sites.{mode}.sorted.col2.trim.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule perl_mode_big_awk:
    input:
        'output/{sample}/{mode}/{sample}.sites.{mode}.sorted.col2.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """


rule seperate_forward_strand:
    input:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.fwd.{mode}.sorted.trim.bed'
    
    shell:'''
    grep "+" {input} > {output}
    '''


rule seperate_reverse_strand:
    input:
        'output/{sample}/{mode}/{sample}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/{mode}/{sample}.rev.{mode}.sorted.trim.bed'
    shell:'''
    grep "-" {input} > {output}
    '''

rule merge_bed_fwd:
    input:
        'output/{sample}/{mode}/{sample}.fwd.{mode}.sorted.trim.bed'
    output:
        temp('output/{sample}/nick_to_nick/{sample}.fwd.{mode}.sorted.trim.merge.bed')
    shell:'''
    bedtools merge -i {input} > {output}
    '''


rule merge_bed_rev:
    input:
        'output/{sample}/{mode}/{sample}.rev.{mode}.sorted.trim.bed'
    output:
        temp('output/{sample}/nick_to_nick/{sample}.rev.{mode}.sorted.trim.merge.bed')
    shell:'''
    bedtools merge -i {input} > {output}
    '''



rule nick_to_nick_distance_fwd:
    input:
        'output/{sample}/nick_to_nick/{sample}.fwd.{mode}.sorted.trim.merge.bed'
    output:
        'output/{sample}/nick_to_nick/{sample}.fwd.{mode}.sorted.trim.bed'
    shell:'''
    bedtools closest -t "first" -io -s -d -D "a" -a {input} -b {input} > {output}
    '''


rule nick_to_nick_distance_rev:
    input:
        'output/{sample}/nick_to_nick/{sample}.rev.{mode}.sorted.trim.merge.bed'
    output:
        'output/{sample}/nick_to_nick/{sample}.rev.{mode}.sorted.trim.bed'
    shell:'''
    bedtools closest -t "first" -io -s -d -D "a" -a {input} -b {input} > {output}
    '''


rule plot_ntn_fwd:
    input:
        'output/{sample}/nick_to_nick/{sample}.fwd.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/nick_to_nick/{sample}.fwd.{mode}.sorted.trim.plt.png'
    shell:'''
    scripts/Rscript {input} {output}
    '''


rule plot_ntn_rev:
    input:
        'output/{sample}/nick_to_nick/{sample}.rev.{mode}.sorted.trim.merge.bed'
    output:
        'output/{sample}/nick_to_nick/{sample}.rev.{mode}.sorted.trim.plt.png'
    shell:'''
    scripts/Rscript {input} {output}
    '''


rule bed_to_bam_fwd:
    input:
        bed='output/{sample}/{mode}/{sample}.fwd.{mode}.sorted.trim.bed',
        index='rawdata/hg19/hg19.chrom.sizes'
    output:
        'output/{sample}/bam/{sample}.fwd.{mode}.sorted.trim.bam'
    params:
        bed_dir='output/{sample}/bam/'
    shell:'''
    mkdir -p {params.bed_dir}
    bedtools bedtobam -i {input.bed} -g "{input.index}" > {output}
    '''


rule bed_to_bam_rev:
    input:
        bed='output/{sample}/{mode}/{sample}.rev.{mode}.sorted.trim.bed',
        index='rawdata/hg19/hg19.chrom.sizes'
    output:
        'output/{sample}/bam/{sample}.rev.{mode}.sorted.trim.bam'
    params:
        bed_dir='output/{sample}/bam/'
    shell:'''
    mkdir -p {params.bed_dir}
    bedtools bedtobam -i {input.bed} -g "{input.index}" > {output}
    '''


rule bam_converage_fwd:
    conda:
        'envs/deeptools.yml'
    input:
        'output/{sample}/bam/{sample}.fwd.{mode}.sorted.trim.bam'
    output:
        'output/{sample}/coverage/{sample}.fwd.{mode}.sorted.trim.bw'
    threads: 12
    params:
        bw_dir='output/{sample}/coverage/'
    shell:'''
    mkdir -p {params.bw_dir}
    samtools index {input}
    bamCoverage --numberOfProcessors {threads} --outFileFormat bigwig --bam {input} -o {output} --samFlagExclude 16
    '''


rule bam_converage_rev:
    conda:
        'envs/deeptools.yml'
    input:
        'output/{sample}/bam/{sample}.rev.{mode}.sorted.trim.bam'
    output:
        'output/{sample}/coverage/{sample}.rev.{mode}.sorted.trim.bw'
    threads: 12
    params:
        bw_dir='output/{sample}/coverage/'
    shell:'''
    mkdir -p {params.bw_dir}
    samtools index {input}
    bamCoverage --numberOfProcessors {threads} --outFileFormat bigwig --bam {input} -o {output} --samFlagInclude 16
    '''



rule closest_drip_breaks_forward_strand:
    input:
        gloe='output/{gloe_sample}/{mode}/{gloe_sample}.fwd.{mode}.sorted.trim.bed',
        drip='output/trim_bed/DRIP/{drip_sample}.trim.merged.bed'
    output:
        'output/drip_overlap/fwd/{drip_sample}_{gloe_sample}.fwd.bed'
    
    shell:'''
    bedtools closest -t "first" -d -D "a" -a {input.drip} -b {input.gloe} > {output}
    '''
    # -iu ignore upstream, -s same strand 


rule closest_drip_breaks_rev_strand:
    input:
        gloe='output/{gloe_sample}/{mode}/{gloe_sample}.rev.{mode}.sorted.trim.bed',
        drip='output/trim_bed/DRIP/{drip_sample}.trim.merged.bed'
    output:
        'output/drip_overlap/rev/{drip_sample}_{gloe_sample}.rev.bed'
    shell:'''
    bedtools closest -t "first" -d -D "a" -a {input.drip} -b {input.gloe} > {output}
    '''


rule rand_closest_drip_breaks_rev_strand:
    input:
        gloe='output/{gloe_sample}/{mode}/{gloe_sample}.rev.{mode}.sorted.trim.bed',
        drip='output/random_bed/{drip_sample}_{gloe_sample}.rev.sorted.bed'
    output:
        'output/drip_overlap/rev/{drip_sample}_{gloe_sample}.rev.rand.closest.bed'
    shell:'''
    bedtools closest -t "first" -d -D "a" -a {input.drip} -b {input.gloe} > {output}
    '''

rule rand_closest_drip_breaks_fwd_strand:
    input:
        gloe='output/{gloe_sample}/{mode}/{gloe_sample}.fwd.{mode}.sorted.trim.bed',
        drip='output/random_bed/{drip_sample}_{gloe_sample}.fwd.sorted.bed'
    output:
        'output/drip_overlap/fwd/{drip_sample}_{gloe_sample}.fwd.rand.closest.bed'
    shell:'''
    bedtools closest -t "first" -d -D "a" -a {input.drip} -b {input.gloe} > {output}
    '''

rule sort_rand_bed_rev:
    input:
       'output/random_bed/{drip_sample}_{gloe_sample}.rev.bed'
    output:
        'output/random_bed/{drip_sample}_{gloe_sample}.rev.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''

rule sort_rand_bed_fwd:
    input:
          'output/random_bed/{drip_sample}_{gloe_sample}.fwd.bed'
    output:
         'output/random_bed/{drip_sample}_{gloe_sample}.fwd.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


# end bam2bedI rules 

rule plot_distance_distrabutions_fwd:
    input:
        bed='output/drip_overlap/fwd/{drip_sample}_{gloe_sample}.fwd.bed',
    output:
        png='output/plots/fwd/{drip_sample}_{gloe_sample}.fwd.png',
        #rds='output/plots/fwd/{drip_sample}_{gloe_sample}.fwd.rds'
    params:
        output_name='output/plots/fwd/{drip_sample}_{gloe_sample}.fwd',
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/plot_nick_distances.R {input} {params.output_name}
    '''

rule plot_distance_distrabutions_rev:
    input:
        bed='output/drip_overlap/rev/{drip_sample}_{gloe_sample}.rev.bed',
    output:
        png='output/plots/rev/{drip_sample}_{gloe_sample}.rev.png',
    params:
        output_name='output/plots/rev/{drip_sample}_{gloe_sample}.rev',
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/plot_nick_distances.R {input.bed} {params.output_name}
    '''

rule plot_rand_distance_distrabutions_rev:
    input:
        rand_bed='output/drip_overlap/rev/{drip_sample}_{gloe_sample}.rev.rand.closest.bed'
    output:
        random_plot='output/plots/rev/{drip_sample}_{gloe_sample}.rev.rand.png'
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/plot_nick_distances.R {input.rand_bed} {output.random_plot}
    '''

rule plot_rand_distance_distrabutions_fwd:
    input:
        rand_bed='output/drip_overlap/fwd/{drip_sample}_{gloe_sample}.fwd.rand.closest.bed'
    output:
        random_plot='output/plots/fwd/{drip_sample}_{gloe_sample}.fwd.rand.png'
        #rds='output/plots/fwd/{drip_sample}_{gloe_sample}.fwd.rds'
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/plot_nick_distances.R {input.rand_bed} {output.random_plot}
    '''

rule concentrations_vs_nick_dist:
    input:
         'output/{sample}/footloop_con/{mode}.footloop_all.{strand}.{sample}.{gloe_strand}.con.sorted.first.closest.bed'
    output:
         'output/{sample}/footloop_con/{mode}.footloop_all.{strand}.{sample}.{gloe_strand}.con.sorted.first.closest.png'
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/concentrations_vs_distance.R {input} {output}
    '''


rule footloop_concentration:
    input:
        amplicons='rawdata/footloop/footprinted_sites.bed',
        footloop='rawdata/footloop/footloop_all.{strand}.bed'
    output:
        'output/footloop_cons/footloop_all.{strand}.con.bed'
    shell:'''
    /home/ethollem/anaconda3/bin/Rscript scripts/footloop_init_concentration.R \
    {input.amplicons} {input.footloop} {output}
    '''


rule sort_footloop_concentrations:
    input:
        'output/footloop_cons/footloop_all.{strand}.con.bed'
    output:
        'output/footloop_cons/footloop_all.{strand}.con.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule closest_con_sites_to_breaks:
    input:
        footloop_cons='output/footloop_cons/footloop_all.{strand}.con.sorted.bed',
        gloe_reads='output/{sample}/{mode}/{sample}.{gloe_strand}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/footloop_con/{mode}.footloop_all.{strand}.{sample}.{gloe_strand}.con.sorted.first.closest.bed'
    params:
        outdir='output/{sample}/footloop_con/'
    shell:'''
    mkdir -p {params.outdir}
    bedtools closest -D ref -t first -a {input.footloop_cons} -b {input.gloe_reads} \
    > {output}
    '''




# This is incorrect for a couple reasons, the first being that we should prob
# only sample footloop sites second has to do with t-tests requiring
# normality use the rule below instead
rule make_random_fl_initiation_sites_bed:
    input:
        'rawdata/footloop/genes.tsv'  # should be produces with R script but having problems getting biomaRt working
    output:
        temp('output/random/fl_init_sites.random.bed')
    shell:'''
    mkdir -p output/random
    /home/ethollem/anaconda3/bin/Rscript scripts/footloop_init_dist_rand.R \
    {input} {output}
    '''


rule random_sample_footloop_sites:
    input:
        'rawdata/footloop/footprinted_sites.bed'
    output:
        temp('output/random/fl.sites.random.bed')
    shell:'''
    mkdir -p output/random
    /home/ethollem/anaconda3/bin/Rscript scripts/random_sample_footloop_sites.R \
    {input} {output}
    '''


rule sort_random_fl_init_sites_bed:
    input:
        'output/random/fl.sites.random.bed'
    output:
        'output/random/fl.sites.random.sorted.bed'
    shell:'''
    bedtools sort -i {input} > {output}
    '''


rule closest_breaks_to_random_init_sites:
    input:
        rand_sites='output/random/fl.sites.random.sorted.bed',
        gloe_reads='output/{sample}/{mode}/{sample}.{strand}.{mode}.sorted.trim.bed'
    output:
        'output/{sample}/random/{mode}.{sample}.{strand}.random.closest.first.bed'
    params:
        outdir='output/{sample}/random'
    shell:'''
    mkdir -p {params.outdir}
    bedtools closest -D ref -t first -a {input.rand_sites} -b {input.gloe_reads} \
    > {output}
    '''


rule plot_footloop_vs_random_distances:
    input:
        random='output/{sample}/random/{mode}.{sample}.{gloe_strand}.random.closest.first.bed',
        footloop='output/truncated_closests/footloop_all.{mode}.footloop_{footloop_strand}.gloe_{gloe_strand}.trunc.closest.first.{sample}.bed'
    output:
        'output/plots/footloop_dist_vs_random/{mode}.{sample}.footloop_{footloop_strand}.gloe_{gloe_strand}.nick_dist.png'
    shell:'''
    mkdir -p output/plots/footloop_dist_vs_random/
    /home/ethollem/anaconda3/bin/Rscript scripts/footloop_compare_rand.R \
    {input.footloop} {input.random} {output}
    '''


rule download_hg18_chr_sizes:
    output:
        'rawdata/hg19/hg19.chrom.sizes'
    shell:'''
    curl -L http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes -o {output}
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


rule download_hg19_tss:
    output:
        'rawdata/hg19/ucsc_promoter_list.sga'
    shell:'''
    curl -L https://ccg.epfl.ch/mga/hg19/ucsc/ucsc_promoter_list.sga.gz -o {output}
    '''


rule download_transcript_ends:
    output:
        'rawdata/hg19/ucsc_transcriptEnd_list.sga'
    shell:'''
    curl -L https://ccg.epfl.ch/mga/hg19/ucsc/ucsc_transcriptEnd_list.sga.gz -o {output}
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

rule merge_DRIP:
    # merge regions of significant signal together to avoid plotting
    # a point for each base of an R loop later on 
    input:
        'output/trim_bed/DRIP/{sample}.trim.bed'
    output:
        'output/trim_bed/DRIP/{sample}.trim.merged.bed'
    params:
        d="50",
    shell:'''
    bedtools merge -d {params.d} -i {input} > {output}
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




