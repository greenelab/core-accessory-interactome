## Run this once to setup environment
## Used R 3.6.3

get_DE_stats_DESeq <- function(metadata_file,
                               expression_file,
						       out_file) {

  # This function performs DE analysis using DESeq.
  # Expression data in expression_file are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # metadata_file: str
  #   File containing mapping between sample id and group
  #
  # expression_file: str
  #   File containing gene expression data
  #
  # out_file: str
  #   File to output DE stats results 
  expression_data <- t(as.matrix(read.csv(expression_data_file, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_file, sep="\t", header=TRUE, row.names=1))

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  #print(head(expression_data))

  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~genotype)
  
  deseq_object <- DESeq(ddset)

  deseq_results <- results(deseq_object)

  deseq_results_df <-  as.data.frame(deseq_results)

  write.table(deseq_results_df, file = out_file, row.names = T, sep = "\t", quote = F)
                              }