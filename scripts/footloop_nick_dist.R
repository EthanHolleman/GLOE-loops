library(tidyr)
library(ggplot2)
library(RColorBrewer)


main <- function(){
    args = commandArgs(trailingOnly=TRUE)
    bed.path <- args[1]
    output.path <- args[2]
    bed.df <- read_bed(bed.path)
    bed.df.dist <- calculate_nick_distances(bed.df)
    plt <- plot_nick_distances_from_rloop_initiation_site(bed.df.dist)
    ggsave(output.path, plt, dpi=500)
}


# read bed file created from intersecting footloop all peaks bed file and
# a corresponding GLOE-seq reads bed file
read_bed <- function(bed.path){
  bed.df <- as.data.frame(read.delim(bed.path, header=FALSE), sep='\t')
  colnames(bed.df) <- c(
    'fl_chr', 'fl_start', 'fl_end', 'fl_read', 'fl_score', 'fl_strand', 
    'gl_chr', 'gl_start', 'gl_end', 'gl_read', 'gl_score', 'gl_strand'
    )
  bed.df <- separate(bed.df, 'fl_read', into = c('fl_gene', 'fl_read'))
  # need to separate col 4 (1 indexed) into gene and molecule
  return(bed.df)
}

# Calculate the distances from GLOE-seq reads to R-loop initiation sites
# add this as a new column called "nick_dist"
calculate_nick_distances <- function(bed.df){
  nick_dists <- list()
  for (i in 1:nrow(bed.df)){
    dist <- bed.df[i, ]$gl_start - bed.df[i, ]$fl_start
    nick_dists[[i]] <- dist
  }
  bed.df.dist <- bed.df
  bed.df.dist$nick_dist <- unlist(nick_dists)
  return(bed.df.dist)
}

plot_nick_distances_from_rloop_initiation_site <- function(bed.df.dist){
  # want to center on zero (R-loop init site) and want to plot the density?
  # Y axis would be number of reads at that location (GLOE-seq reads)
  # limit to specific range 
  # although is that is the case probably should have done closest not inter
  # section since we still want to consider GLOE-seq reads outside of the R-looops
  # but data will end up looking the same so just use this data to get the
  # plots working
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(bed.df.dist$fl_gene)))
  # combined_hist <- ggplot(bed.df.dist, aes(x=nick_dist, fill=fl_gene)) + 
  #   geom_histogram(binwidth = 20, alpha=0.7, color="black") + 
  #   scale_fill_manual(values = mycolors) + theme_minimal() + 
  #   labs(y='GLOE-seq read count', x='Distance from R-loop initiation site')
  facet_hist <- ggplot(bed.df.dist, aes(x=nick_dist, fill=fl_gene)) + 
    geom_histogram(binwidth = 20, alpha=0.9) + 
    scale_fill_manual(values = mycolors) + theme_minimal() + 
    labs(y='GLOE-seq read count', x='Distance from R-loop initiation site') +
    facet_wrap(~fl_gene) + theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(trans='log10')
  
  facet_hist


}


if (!interactive()){
  main()
}
