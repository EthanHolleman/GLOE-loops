# Read in bedfile of footloop data and convert the end position
# of the peak for each read to 1 bp downstream of the initiation site.
# This can then be used to find distances of other genomic features
# (nicks) to R-loop initiation sites.


read_bed_file <- function(bed.path){
  
  bed.df <- as.data.frame(read.delim(bed.path, header=FALSE), sep='\t')
  return(bed.df)

}


modify_peak_end_position <- function(bed.df){
  
  bed.df$V3 <- bed.df$V2 + 1
  stopifnot(all(bed.df$V2 - bed.df$V3 == -1))
  return(bed.df)
  
}


put_df_to_bed <- function(bed.df, output.path){
  
  write.table(bed.df, output.path, sep = '\t',
              row.names=F, quote = F, col.names = F
  )
  return(output.path)
  
}


main <- function(){
  
  args <- commandArgs(trailingOnly=TRUE)
  input_bed <- args[1]
  output_bed <- args[2]
  bed.df <- read_bed_file(input_bed)
  bed.df.trunk <- modify_peak_end_position(bed.df)
  put_df_to_bed(bed.df.trunk, output_bed)
  
}


if (!interactive()){
  main()
}