# Read the foot loop all peaks file and create a tsv file of
# gene coordinates
library(biomaRt)
library(tidyr)


get_hg_mart <- function(){

  ensembl.human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                           host="feb2014.archive.ensembl.org", 
                           path="/biomart/martservice", 
                           dataset="hsapiens_gene_ensembl")
  return(ensembl.human)
}

get_hg_genes <- function(mart, gene_names, 
                         attributes= c('hgnc_symbol', 
                                       'chromosome_name', 
                                       'start_position', 
                                       'end_position')){
  
  genes <- getBM(attributes = attributes, filters = 'hgnc_symbol',
                 values=gene_names, mart = mart)
  return(genes)
  
}


write_genes <- function(genes.df, output.path){
  
  write.table(genes.df, output.path, sep = '\t',
              row.names=F, quote = F, col.names = T
  )
  return(output.path)
}


cononical_chrs <- function(genes.df){
  
  cc <- c(c(1:22), "X", "Y")
  return(subset(genes.df, chromosome_name %in% cc))
  
}

read_footloop_bed <- function(footloop.path){
  
  bed.df <- as.data.frame(read.delim(footloop.path, header=FALSE), sep='\t')
  colnames(bed.df) <- c('fl_chr', 'fl_start', 'fl_end', 'fl_read', 'fl_score', 
                        'fl_strand')
  bed.df <- separate(bed.df, 'fl_read', into = c('fl_gene', 'fl_read'))
  return(bed.df)
  
}

main <- function(){
  
  args = commandArgs(trailingOnly=TRUE)
  footloop.path <- args[1]
  output.path <- args[2]
  # footloop.path <- '/home/ethan/data/igg/chedin/footloop/footloop_peak_all.bed'
  # output.path <- '/home/ethan/data/igg/chedin/footloop/genes.tsv'
  footloop.df <- read_footloop_bed(footloop.path)
  message('Getting Mart')
  mart <- get_hg_mart()
  message('Getting Genes')
  genes.df <- get_hg_genes(mart, unique(footloop.df$fl_gene))
  genes.df <- cononical_chrs(genes.df)
  message('Writing Genes')
  write_genes(genes.df, output.path)
  genes.df
}

if (! interactive()){
  main()
}