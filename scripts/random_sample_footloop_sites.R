# Pass in a bed file containing footloop coordinates generated locations within the genes specified in that file.
# This can then be used for comparing R-loop initiation site distances
# to random locations in the genome.

# Then get the closest break using bedtools
# Then process this file and do random sampeling to create normal distrabutions
# for both variables being compared and finally plot



read_bed <- function(bed.path){
  
  return(as.data.frame(read.delim(bed.path, header = F)))
  
}


random_gene_init_sites <- function(gene_row, total_samples=5000){
  
  rand_init_sites_start <- sample(gene_row$V2:gene_row$V3,
                                  total_samples, replace = T)
  rand_init_sites_end <- rand_init_sites_start + 1
  
  cbind(
    rep(gene_row$V1, total_samples),
    rand_init_sites_start,
    rand_init_sites_end,
    paste(gene_row$V4, 1:total_samples, sep='_')
  )
  
  
}


random_sites_all_genes <- function(bed.df, total_samples_per_gene=5000){
  
  random_init_sites.list <- list()
  for (i in 1:nrow(bed.df)){
    random_init_sites.list[[i]] <- random_gene_init_sites(
      bed.df[i, ], total_samples_per_gene
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
  bed.path <- args[1]
  output.path <- args[2]
  bed.df <- read_bed(bed.path)
  random_init_sites.df <- random_sites_all_genes(bed.df)
  write_random_sites_as_bed(random_init_sites.df, output.path)
  
}


if (! interactive()){
  main()
}




