####################### DESEQ2 COUNT NORMALISATION FROM HTSEQ-COUNT OUTPUT ###########################

# setRepositories(ind=1:10)
# install.packages("DESeq2",dependencies = T)
# install.packages("foreign",dependencies = T)
# install.packages("optparse",dependencies=T)
library(DESeq2)

DESeq2prepare = function(htseqdir, pattern, metacols = 1, sample_info_file, expt_label){

# FIND HTSEQ-COUNT OUTPUT FILES

htseqcount_files <- dir(htseqdir, pattern = pattern,full.names = T)
print(htseqcount_files)

# LOCATE HTSEQ-COUNT SUMMARY STATISTICS ROWNAMES IN DATA FOR FUTURE SEPARATION

htseq_rows <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")

# FOR EACH HTSEQ-COUNT FILE, READ COUNTS
print("READING COUNTS")

lapply(htseqcount_files, function(x){
  
  counts <- read.delim(x, header = F,sep = "\t", stringsAsFactors = F, col.names = c("Feature",paste0(x,".Count")))
  
  return(counts)
}) -> count_list
names(count_list) = htseqcount_files

str(count_list)
#count_matrix = do.call(cbind, count_list)
count_matrix = Reduce(f = function(x,y)merge(x, y, by = "Feature"), count_list)
str(count_matrix)

# WRITE RAW COUNT FILES, SUMMARY STATISTICS

write.table(count_matrix, file=paste(expt_label,"raw.counts.txt",sep=""), sep="\t", quote=F, col.names = NA)

# LOAD SAMPLE INFORMATION FILE

sample_info <- read.delim(sample_info_file, header = T,stringsAsFactors = F,row.names=1)
sample_info <- sample_info[colnames(count_matrix),]

# LOADING DATA INTO DESEQ2 COUNT MATRIX OBJECT

rownames(count_matrix) <- annot
deseq2_matrix <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info,design = ~1)

# RLOG TRANSFORMATION (CORRECTION FOR INFLATION OF FOLD CHANGES AT SMALL COUNT VALUES, NORMALISATION WITH RESPECT TO LIBRARY SIZE)

assay(rlog(deseq2_matrix)) -> rlog_norm_counts
cbind(annot[rownames(rlog_norm_counts),],rlog_norm_counts) -> rlog_norm_counts_with_gene_annotation

filler <- matrix("",nrow = ncol(sample_info),ncol = ncol(annot),dimnames = list(colnames(sample_info),c("Name","Symbol","Description")))
rlog_file_header <- cbind(filler, t(sample_info))

#rlog_colnames <- c("Name","Symbol","Description",colnames(rlog_norm_counts))
rlog_norm_counts_with_all_annotation <- t(data.frame(t(rlog_file_header), t(rlog_norm_counts_with_gene_annotation),stringsAsFactors = F))

# WRITE RLOG NORMALISED DATA TO FILE

write.table(rlog_norm_counts_with_all_annotation, file=paste("rLog.normalised.counts.",opt$outfile,sep=""), sep="\t", quote=F, col.names = NA)

}


DESeq2NormFromHTSeqCount = function(directory, pattern, design, sample_groups = "", expt_label){
  
  htseqcount_files <- dir(directory, pattern = pattern,full.names = F)
  print(htseqcount_files)
  
  sampleTable = data.frame(sample_name = htseqcount_files, file_name = htseqcount_files, class = sample_groups)
  print(sampleTable)
  
  deseq2_matrix <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~1)
  str(deseq2_matrix)
  assay(rlog(deseq2_matrix)) -> rlog_norm_counts
  
  write.table(rlog_norm_counts, file=paste(expt_label, ".rLog.normalised.counts.txt",sep=""), sep="\t", quote=F, col.names = NA)
  
}

#DESeq2NormFromHTSeqCount(directory = "../RNAseq_counts/", pattern = "COUNT", design = ~1, expt_label = "AS_RNAseq")
