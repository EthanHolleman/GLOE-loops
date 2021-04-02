library(ggplot2)
library(ggpubr)
library(tidyr)


read_bed <- function(bed.path){
  
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

# relationship score and closest nick

plot_concentration_vs_nick_dist <- function(bed.df){
  
  bed.df <- subset(bed.df, fl_score > 1)
  plt.line <- ggplot(bed.df, aes(y=fl_score, x=abs(nick_dist), color=fl_gene)) + 
    theme_pubclean() + 
    labs(x='Distance to closest nick', y='Proximal initiation sites', color='Gene') +
    geom_smooth(method=lm)
  plt.dots <-  ggplot(bed.df, aes(y=fl_score, x=abs(nick_dist), color=fl_gene), alpha=0.5) + 
    geom_point() + theme_pubclean() + 
    labs(x='Distance to closest nick', y='Proximal initiation sites', color='Gene')
  
  ggarrange(plt.line, plt.dots, common.legend = TRUE, legend = 'bottom')
}


save_plot <- function(plt, output.path){
  
  ggsave(output.path, plt, dpi=500)
  
}



main <- function(){
  
  args = commandArgs(trailingOnly=TRUE)
  footloop.con.closest.bed <- args[1]
  output.path <- args[2]
  footloop.con.closest.bed.df <- read_bed(footloop.con.closest.bed)
  plt <- plot_concentration_vs_nick_dist(footloop.con.closest.bed.df)
  save_plot(plt, output.path)
  
}

if (! interactive()){
  main()
}



