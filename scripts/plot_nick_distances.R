library(ggplot2)

main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  print(args)
  distances <- process_bed(args[1])
  d_plot <- plot_distances(distances)
  save_plot(args[2], d_plot)
}


process_bed <- function(bed.path){
  
  con <- file(description=bed.path, open="r") 
  distances <- list()
  index <- 1
  while (TRUE){
    line <- read.delim(con, nrows = 1, header = F)
    distances[[index]] <- line[length(line)]
    index <- index + 1
  }
  return(distances)
}


plot_distances <- function(distances){
  
  dist.df <- as.data.frame(unlist(distances))
  colnames(dist.df) <- 'distance'
  d_plot <- ggplot(dist.df, aes(x=distance)) + 
    geom_histogram(color='black', fill='royalblue', alpha=0.4) +
    labs(x='Distance to closest nick', y='R-loop count') + 
    theme_minimal()
  return(d_plot)
}


save_plot <- function(output_path, d_plot){
  op_no_suffix <- tools::file_path_sans_ext(output_path)
  png_path <- paste(op_no_suffix, '.png', sep = '')
  rds_path <- paste(op_no_suffix, '.rds', sep='')
  ggsave(png_path, d_plot, dpi=500)
  saveRDS(d_plot, rds_path)
}

if (!interactive()){
  main()
}