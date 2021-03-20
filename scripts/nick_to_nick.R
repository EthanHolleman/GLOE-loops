# Read a bed file that compares the distances between nicks
# on a given strand (using bedtools closest and ignoring
# distances of nick A -> nick A)
library(ggplot2)
library(ggpubr)


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  input_bed <- args[1]
  output_path <- args[2]
  bed.df <- process_bed(input_bed)
  dist.plt <- plot_distances(bed.df)
  ggsave(output_path, dist.plt, dpi=500)
}


process_bed <- function(bed.path){
  con <- file(description=bed.path, open="r") 
  distances <- list()
  chromosomes <- list()
  index <- 1
  while (TRUE){
    line <- tryCatch(read.delim(con, nrows = 1, header = F), error=function(e) NULL)
    if (!is.null(line)){
      distances[[index]] <- line[length(line)]
      chromosomes[[index]] <- line[[1]]
      index <- index + 1
    }
    else{
      break
    }
    
  }
  close(con)
  bed.df <- do.call(rbind, Map(data.frame, A=distances, B=chromosomes))
  colnames(bed.df) <- c('distance', 'chromosome')
  return(bed.df)
}



plot_distances <- function(bed.df){
  
  dist_plot <- ggplot(bed.df, aes(x=distance), show.legend = FALSE) + 
    geom_density(fill='royalblue', alpha=0.7) + theme_minimal() + labs(x='')
  dist_box <- ggplot(bed.df, aes(x=distance)) + geom_boxplot(fill='royalblue') + 
    theme_minimal()
  combine <- ggarrange(dist_plot, dist_box, 
            ncol = 1, nrow = 2,  align = "hv", 
            widths = c(1, 1), heights = c(2, 1))
  return(combine)
}


if (!interactive()){
  main()
}else{
  bed.path <- '/home/ethan/data/igg/chedin/HeLa_DRIPc_WT_LS61A_rep1_pos.bw_GSM3939125.fwd.bed'
  bed.df <- process_bed(bed.path)
  combine <- plot_distances(bed.df)
  ggsave('test.png', combine, dpi=500)
  
}


