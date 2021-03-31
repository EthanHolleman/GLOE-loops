library(biomaRt)


get_hs_mart <- function(){
  ensembl <- useMart("ensembl")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  return(ensembl.human)
}

get_hg_genes <- function(mart, gene_names, 
                        attributes= c('hgnc_symbol', 
                                      'chromosome_name', 
                                      'start_position', 
                                      'end_position')){
  
  genes <- getBM(attributes = attributes, filters = 'hgnc_symbol',
                 values=gene_names, mart = mart)
  return(genes)

}

# then need to get length of genes and the number of things in
# them (breaks) locations then randomly place the r-loop peaks in the
# the gene and get the distance to the nearest peak
# create a distrabution for that
# There are basically two questions are we looking for R-loop
# enrichment or breaks nearby?