library(ggpubr)
library(ggplot2)
library(tidyr)



read_footloop_bed <- function(bed.path){
  bed.df <- as.data.frame(read.delim(bed.path, header=FALSE), sep='\t')
  colnames(bed.df) <- c(
    'fl_chr', 'fl_start', 'fl_end', 'fl_read', 'fl_score', 'fl_strand', 
    'gl_chr', 'gl_start', 'gl_end', 'gl_read', 'gl_score', 'gl_strand',
    'nick_dist'
  )
  bed.df <- separate(bed.df, 'fl_read', into = c('fl_gene', 'fl_read'))
  # need to separate col 4 (1 indexed) into gene and molecule
  return(bed.df)
}

read_random_bed <- function(bed.path){
  
  bed.df <- as.data.frame(read.delim(bed.path, header=FALSE), sep='\t')
  colnames(bed.df) <- c(
    'fl_chr', 'fl_start', 'fl_end', 'fl_read', 
    'gl_chr', 'gl_start', 'gl_end', 'gl_read', 'gl_score', 'gl_strand',
    'nick_dist'
  )
  bed.df <- separate(bed.df, 'fl_read', into = c('fl_gene', 'fl_read'))
  return(bed.df)
  
}


get_distance_info_from_bed <- function(bed.df){
  # only get information for comparing distance
  return(bed.df[, c(1:4, ncol(bed.df))])
}


label_df <- function(bed.df, label){
  
  bed.df$group <- rep(label, nrow(bed.df))
  return(bed.df)
  
}

process_king_bed <- function(king.bed){
  
  genes <- unique(king.bed$fl_gene)
  
  for (i in 1:length(genes)){
    
    footloop_reads <- nrow(subset(king.bed, fl_gene==genes[i] & group=='footloop'))
    rand_reads <- nrow(subset(king.bed, fl_gene==genes[i] & group=='random'))

    if (footloop_reads <= 10 | rand_reads <= 10){
      king.bed <- king.bed[!(king.bed$fl_gene)==genes[i], ]
    }
    
  }
  
  king.bed$nick_dist <- abs(king.bed$nick_dist)
  return(king.bed)
   
}

plot_compare_means <- function(king.bed){
  
  plt <- ggplot(king.bed, aes(x=group, y=nick_dist, fill=group)) + 
    geom_boxplot() + 
    stat_compare_means(method='t.test', label.y = 6000) + 
    scale_fill_manual(values=c('royalblue', 'tan')) + 
    facet_wrap(~fl_gene) + 
    labs(x='', y='Distance to Closest Nick') + 
    theme_pubclean()
  
  return(plt)
  
}

save_plot <- function(plt, output.path){
  
  ggsave(output.path, plt, dpi=1000, limitsize = F, height = 10, 
         width=10, units = 'in')
  
}


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  footloop.bed <- args[1]
  random.bed <- args[2]
  output.path <- args[3]
  footloop.bed.df <- read_footloop_bed(footloop.bed)
  random.bed.df <- read_random_bed(random.bed)
  footloop.bed.df <- get_distance_info_from_bed(footloop.bed.df)
  random.bed.df <- get_distance_info_from_bed(random.bed.df)
  footloop.bed.df <- label_df(footloop.bed.df, 'footloop')
  random.bed.df <- label_df(random.bed.df, 'random')
  king.bed <- rbind(footloop.bed.df, random.bed.df)
  king.bed <- process_king_bed(king.bed)
  plt <- plot_compare_means(king.bed)
  save_plot(plt, output.path)
}

if (! interactive()){
  main()
}



