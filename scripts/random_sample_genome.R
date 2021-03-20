

main <- function(){
  options(scipen=100000)  # disable sci notation
  args = commandArgs(trailingOnly=TRUE)
  chrom.sizes.path <- args[1]
  output.path <- args[2]
  n_samples <- args[3]
  
  chr.df <- read_chr_sizes(chrom.sizes.path)
  chr.df.canonical <- get_canonical_chr(chr.df)
  random_loci.df <- sample_genome(chr.df.canonical, n_samples)
  write_rand_loci_as_bed(random_loci.df, output.path)
  
}

read_chr_sizes <- function(chr.sizes.path){
  
  chr.df <- as.data.frame(read.delim(chr.sizes.path, header=F))
  colnames(chr.df) <- c('chromosome', 'size')
  return(chr.df)
}

get_canonical_chr <- function(chr.df){
  # get canonical chromosomes by filtering by size
  subset(chr.df, size >= 48129895)  # size of 21

}

random_genomic_locus <- function(chr.df){
  # pick a random 1 base pair locus in the genome and return as
  # vector which can be written as bed
  chr.row <- sample(1:nrow(chr.df), 1)
  chr <- chr.df[chr.row, 1]
  max_position <- chr.df[chr.row, ]$size
  locus_start <- sample(1:max_position, 1)
  locus_end <- locus_start + 1
  return(c(chr, locus_start, locus_end))
}

sample_genome <- function(chr.df, n){
  loci <- list()
  for (i in 1:n){
    loci[[i]] <- random_genomic_locus(chr.df)
  }
  rand_loci.df <- as.data.frame(do.call(rbind, loci))
  print(head(rand_loci.df))
  return(rand_loci.df)
}


write_rand_loci_as_bed <- function(rand_loci.df, output.path){
  
  write.table(rand_loci.df, output.path, sep = '\t',
              row.names=F, quote = F, col.names = F
              )
  return(output.path)
}

if (!interactive()){
  main()
}


