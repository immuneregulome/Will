##how to make a counts file from DESeq2
library(DESeq2)
library(writexl)

#load the data in

directory <- "/Users/montgomerywf/Desktop/Sen Colab Files/NKT_RNA_seq/RNAseq_2018/reverse_counts/"

#create some variables to make loading in the dds object prettier
samplefiles <- list.files(directory)
sampleCondition <- str_remove(samplefiles, "_s[0-9].counts_r")
sampletable <- data.frame(sampleName = samplefiles,
                          fileName = samplefiles,
                          condition = sampleCondition)

#create the dds object that we'll use for all downstream analysis
dds_2018 <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                       directory = directory,
                                       design = ~ condition)

#prefiltering
keep <- rowSums(counts(dds_2018)) >= 10
dds_2018 <- dds_2018[keep, ]

#normalize the data
dds_2018 <- DESeq(dds_2018)

#export the counts to a matrix object
counts_df <- counts(dds_2018, normalized = TRUE)
counts_df <- as.data.frame(counts_df, row.names = rownames(counts_df))
counts_df$Genes <- rownames(counts_df)
write_xlsx(x = counts_df, path = "~/Desktop/Sen Colab Files/NKT_RNA_seq/RNAseq_2018/reverse_results/RNAseq_2018_reverse_counts.xlsx",
           format_headers = TRUE)
