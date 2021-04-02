library(ggpubr)
library(ggplot2)
library(tidyr)


# need to add sampling to this script basically

t <- '/home/ethan/data/direct.GSM3939125.fwd.random.closest.first.bed'
t2 <- '/home/ethan/data/footloop_all.direct.footloop_fwd.gloe_fwd.trunc.closest.GSM3939125.bed'
t2 <- '/home/ethan/data/footloop_all.direct.footloop_fwd.gloe_rev.trunc.closest.first.GSM3939125.bed'

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





random_sample_by_gene <- function(bed.df, group_name, n_observations=10, n_samples=500){
  # randomly sample distances to create a distribution
  
  mean_dists <- list()
  genes <- unique(bed.df$fl_gene)
  k = 1
  for (g in 1:length(genes)){
    
    gene_name <- genes[[g]]
    bed.df.gene <- subset(bed.df, fl_gene==gene_name)
    n_rows <- nrow(bed.df.gene)
    for (s in 1:n_samples){
      
      sample_rows <- bed.df.gene[sample(1:n_rows, n_observations), ]
      mean_dist <- mean(sample_rows$nick_dist)
      mean_dists[[k]] <- c(mean_dist, gene_name, group_name)
      k <- k + 1
      
    }
    
  }
  
  dist.df <- data.frame(do.call(rbind, mean_dists))
  colnames(dist.df) <- c('sample_nick_dist', 'fl_gene', 'group')
  dist.df$sample_nick_dist <- as.numeric(dist.df$sample_nick_dist)
  dist.df
  
}


compare_random_samples <- function(bed.df.sample.random, bed.df.sample.footloop){
  
  big.df <- rbind(bed.df.sample.random, bed.df.sample.footloop)
  
  
  
  
  
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
    scale_fill_manual(values=c('royalblue', 'tan')) + 
    labs(x='', y='Distance to Closest Nick') 
  
  y_label_height <-  layer_scales(plt)$y$range$range[[2]] * 0.75
  plt <- plt + stat_compare_means(method='t.test', label.y = y_label_height) + 
        facet_wrap(~fl_gene) + theme_pubclean()
  
  return(plt)
  
}

plot_overall_dists <- function(king.bed){
  
  footloop_mean <- mean(subset(king.bed, group=='footloop')$nick_dist)
  random_mean <- mean(subset(king.bed, group=='random')$nick_dist)
  xlim = footloop_mean * 4
  
  
  plt <- ggplot(king.bed, aes(x=nick_dist, fill=group)) + 
    geom_density(alpha=0.5) + xlim(0, 2000) + 
    scale_fill_manual(values=c('royalblue', 'tan')) + 
    theme_pubclean()
  
  max_y <- layer_scales(plt)$y$range$range[[2]]
  print(max_y)
  
  plt <- plt + geom_vline(xintercept = footloop_mean) + 
         annotate("text", x=footloop_mean + 200, 
             y=max_y * 0.750, 
             label=paste('Footloop mean', round(footloop_mean, 1))) +
        geom_vline(xintercept = random_mean) + 
        annotate("text", x=random_mean + 500, 
               y=max_y * 0.90, 
               label=paste('Random mean', round(random_mean, 1))) +
    labs(y='Density', x='Distance to Closest Nick')
      
  plt
}

save_plot <- function(plt, output.path){
  
  ggsave(output.path, plt, dpi=1000, limitsize = F, height = 12, 
         width=10, units = 'in')
  
}

qqplot <- function(data, color){
  
  qq <- ggqqplot(data, color=color)
  print(length(data))
  label_height <- layer_scales(qq)$y$range$range[[2]] * 0.75

  sample_size <- length(data)
  if (sample_size > 5000){
    sample_size <- 5000
  }

  sample_rows <- sample(1:length(data), sample_size)
  s.t <- shapiro.test(data[sample_rows])
  qq <- qq + annotate('text', x=0, y=label_height, label=paste('p =', s.t$p.value))
  qq
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
  plt1 <- plot_compare_means(king.bed)
  plt2 <- plot_overall_dists(king.bed)
  qq1 <-  qqplot(subset(king.bed, group=='footloop')$nick_dist, color='royalblue')
  qq2 <-  qqplot(subset(king.bed, group=='random')$nick_dist, color='tan')
  big_plt <- ggarrange(
    plt1, plt2, 
    ggarrange(qq1, qq2, nrow=1, ncol=2, heights=c(0.5, 0.5), labels='C'),
    nrow = 3, ncol=1, heights = c(2, 0.7, 0.7), labels=c('A', 'B'))
  save_plot(big_plt, output.path)
}

if (! interactive()){
  main()
}



