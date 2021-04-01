# Pass in a tsv file containing gene coordinates generated from biomaRt
# and randomly sample 1bp locations within the genes specified in that file.
# This can then be used for comparing R-loop initiation site distances
# to random locations in the genome.


read_genes <- function(genes.tsv){
  
  return(as.data.frame(read.delim(genes.tsv)))
  
}

random_gene_init_sites <- function(gene_row, samples=2500){
  
  rand_init_sites_start <- sample(gene_row$start_position:gene_row$end_position,
                            samples, replace = T)
  rand_init_sites_end <- rand_init_sites_start + 1
  
  cbind(
    rep(paste('chr', gene_row$chromosome_name, sep=''), samples),
    rand_init_sites_start,
    rand_init_sites_end,
    paste(gene_row$hgnc_symbol, 1:samples, sep=',')
  )
  
  
}


random_sites_all_genes <- function(genes.df, samples_per_gene=2500){
  
  random_init_sites.list <- list()
  for (i in 1:nrow(genes.df)){
    random_init_sites.list[[i]] <- random_gene_init_sites(
      genes.df[i, ], samples_per_gene
    )
  }
  random_init_sites.df <- as.data.frame(do.call('rbind', random_init_sites.list))
  
}


write_random_sites_as_bed <- function(random_init_sites.df, output.path){
  
  # reorder the columns into bed format
  #random_init_sites.df <- random_init_sites.df[, c(2:4, 1)]
  write.table(random_init_sites.df, output.path, sep = '\t',
              row.names=F, quote = F, col.names = F
  )
  return(output.path)
  
  
}

main <- function(){
  
  args <- commandArgs(trailingOnly=TRUE)
  genes.tsv <- args[1]
  output.path <- args[2]
  genes.df <- read_genes(genes.tsv)
  random_init_sites.df <- random_sites_all_genes(genes.df)
  write_random_sites_as_bed(random_init_sites.df, output.path)

}

if (! interactive()){
  main()
}


# genes <- '/home/ethan/data/igg/chedin/footloop/genes.tsv'
# genes.df <- read_genes(genes)
# l <- random_sites_all_genes(genes.df)
# write_random_sites_as_bed(l, 'test.1.bed')



