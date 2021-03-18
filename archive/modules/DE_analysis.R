## Run this once to setup environment
## Used R 3.6.3
library("DESeq2")

get_DE_stats_DESeq <- function(metadata_filename,
                               expression_filename,
                               out_filename) {

  # This function performs DE analysis using DESeq.
  # Expression data in expression_file are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # metadata_filename: str
  #   File containing mapping between sample id and group
  #
  # expression_filename: str
  #   File containing gene expression data
  #
  # out_filename: str
  #   File containing DE stats

  expression_data <- t(as.matrix(read.csv(expression_filename, sep="\t", header=TRUE, row.names=1)))
  metadata <- as.matrix(read.csv(metadata_filename, sep="\t", header=TRUE, row.names=1))

  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), rownames(metadata)))

  group <- interaction(metadata[,1])

  mm <- model.matrix(~0 + group)

  #print(head(expression_data))

  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)

  deseq_object <- DESeq(ddset)

  # Note parameter settings:
  # `independentFilter=False`: We have turned off the automatic filtering, which
  # filter filter out those tests from the procedure that have no, or little
  # chance of showing significant evidence, without even looking at their test statistic.
  # Typically, this results in increased detection power at the same experiment-wide
  # type I error, as measured in terms of the false discovery rate.
  # cooksCutoff=True (default): Cook's distance as a diagnostic to tell if a single sample
  # has a count which has a disproportionate impact on the log fold change and p-values.
  # These genes are flagged with an NA in the pvalue and padj columns
  deseq_results <- results(deseq_object, independentFiltering=FALSE)

  deseq_results_df <-  as.data.frame(deseq_results)

  write.table(deseq_results_df, file = out_filename, row.names = T, sep = "\t", quote = F)
}
