library(tidyr)
library(ggplot2)
library(ggpubr)

read_bed <- function(bed.path){
  
  # read footloop sites into dataframe
  bed.df <- as.data.frame(read.delim(bed.path, header=FALSE), sep='\t')
  colnames(bed.df) <- c(
    'fl_chr', 'fl_start', 'fl_end', 'fl_read', 'fl_score', 'fl_strand'
    )
  bed.df <- separate(bed.df, 'fl_read', into = c('fl_gene', 'fl_read'))
  return(bed.df)
  
}


read_amplicon_bed <- function(amplicon.bed){
  
  return(as.data.frame(read.delim(bed.path, header = F)))
  
}


calculate_init_density <- function(footloop.df.gene, position, reach=3){
  
  count <- 0
  for (i in -reach:reach){
    reach_pos <- position + i
    count <- count + nrow(subset(footloop.df.gene, fl_start==reach_pos))
  }
  total_footloops <- nrow(footloop.df.gene)
  return(count)
  
}

# distance_to_closest_init <- function(footloop.df.gene, position){
#   
#   dist <- 0
#   while (TRUE){
#     
#     # have to account for distance 
#     
#     
#     
#   }
#   
#   
# }

init_concentration <- function(footloop.df, amplicon.df, reach=3){
  
  # iterate through each amplicon
  concentrations <- c()
  
  for (a in 1:nrow(amplicon.df)){
    
    gene <- amplicon.df[a, ]$fl_gene
    start <- amplicon.df[a, ]$fl_start
    end <- amplicon.df[a, ]$fl_end
    length <- end - start
    message(gene)
    footloop.df.gene <- subset(footloop.df, fl_gene==fl_gene)
    gene_densities <- c()
    
    for (i in 1:length){
      position <- start + i
      gene_densities <- c(
        gene_densities, 
        calculate_init_density(footloop.df.gene, position, reach)
      )
    }
    index <- 1:length(gene_densities)
    position <- start:(end-1)
    gene <- rep(gene, length(gene_densities))
    total_footloops <- rep(nrow(footloop.df.gene), length(gene_densities))
    chr <- rep(footloop.df.gene$fl_chr[[1]], length(gene_densities))
    strand <- rep(footloop.df.gene$fl_strand[[1]], length(gene_densities))
    
    
    densities.df <- data.frame(
      gene_densities,
      index,
      gene,
      total_footloops,
      position,
      chr,
      strand

    )
    concentrations[[a]] <- densities.df
  }
  
  do.call('rbind', concentrations)
  
}

write_con_bed <- function(concentrations.df, output.path){
  
  # write bed file of R-loop init concentrations
  bed.df <- as.data.frame(cbind(
    concentrations.df$chr,
    concentrations.df$position,
    concentrations.df$position + 1,
    paste(concentrations.df$gene, 1:length(concentrations.df$gene), sep='_'),
    concentrations.df$gene_densities,
    concentrations.df$strand
    )
  )
  print(colnames(bed.df))
  bed.df <- subset(bed.df, V5 != 0)
  write.table(bed.df, output.path, row.names = F, col.names = F, quote = F, sep='\t')
  return(output.path)
}

save_plot <- function(plt, output.path){
  
  ggsave(output.path, plt, dpi=500)
  
}


plot_concentration <- function(cons.df){
  
  plt1 <- ggplot(cons.df, aes(x=index, y=gene_densities, color=gene), alpha=0.97) + 
    geom_line() + theme_pubclean() + 
    labs(y='R-loop Initiation Count', x='Position in Amplicon') + 
    theme(legend.position="bottom")
  plt1
  
}

main <- function(){
  
  args = commandArgs(trailingOnly=TRUE)
  amplicons.bed <- args[1]
  footloop.bed <- args[2]
  output.path <- args[3]
  
  output.path.png <- paste(output.path, '.png', sep='')
  amplicons.df <- read_bed(amplicons.bed)
  footloop.df <- read_bed(footloop.bed)
  cons.df <- init_concentration(footloop.df, amplicons.df)
  plt <- plot_concentration(cons.df)
  save_plot(plt, output.path.png)
  write_con_bed(cons.df, output.path)
  
}

if (! interactive()){
  
  main()
  
}


