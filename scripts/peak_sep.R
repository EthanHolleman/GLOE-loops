# input a bedfile (bed graph) and and output path and return a bedgraph file
# that is pruned based on scores. Regions who score below 3 standard
# deviations above the mean are excluded

read_bed_file <- function(bedfile.path){
  bed <- as.data.frame(read.table(bedfile.path), sep='\t', header = F)
  colnames(bed) <- c('chr', 'start', 'stop', 'signal')
  return(bed)
}

bed_mean_signal <- function(bed.df){
  mean(bed.df$signal)
}

seperate_peaks <- function(bed.df){
  mean <- mean(bed.df$signal)
  sd <- sd(bed.df$signal) * 3
  min_signal_for_call <- mean + sd
  subset(bed.df, signal >= min_signal_for_call)
}

write_bed_file <- function(bed.df, output.path){
  write.table(bed.df, file=output.path,
              row.names=FALSE, col.names=FALSE, sep='\t')
}

if (!interactive()){
  
  args = commandArgs(trailingOnly=TRUE)
  bedfile.path <- args[1]
  output.path <- args[2]
  bedfile.df <- read_bed_file(bedfile.path)
  bedfile.df.sub <- seperate_peaks(bedfile.df)
  write_bed_file(bedfile.df.sub, output.path)
}