####################### DESEQ2 COUNT NORMALISATION FROM HTSEQ-COUNT OUTPUT ###########################

# setRepositories(ind=1:10)
# install.packages("DESeq2",dependencies = T)
# install.packages("foreign",dependencies = T)
# install.packages("optparse",dependencies=T)
library(DESeq2)
library(optparse)


################## PARSE COMMAND LINE ARGUMENTS ######################

option_list = list(
  make_option(c("-d", "--htseqdir"), type="character", default=NULL, 
              help="directory containing htseq-count files", metavar="character"),
  make_option(c("-m", "--metacols"), type="numeric", default=NULL, 
              help="number of metadata columns for features", metavar="character"),			  
  make_option(c("-p", "--pattern"), type="character", default="", 
              help="pattern to identify htseq-count files", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="count_matrix.txt", 
              help="filename to write the count matrix", metavar="character"),
  make_option(c("-s", "--sampleinfofile"), type="character", default=NULL, 
              help="file containing sample information", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$htseqdir)){
  stop("At least one argument must be supplied (input directory).\n", call.=FALSE)
}

# FIND HTSEQ-COUNT OUTPUT FILES

htseqcount_files <- dir(opt$htseqdir, pattern = opt$pattern,full.names = T)
print(htseqcount_files)

# LOAD ONE OF THE HTSEQ-COUNT FILES AND SAVE THE GENE ANNOTATION

annot <- read.delim(htseqcount_files[1], header = F,sep = "\t", stringsAsFactors = F,row.names = 1)[,1:3]
colnames(annot) <- c("Name","Symbol","Description")

# LOCATE HTSEQ-COUNT SUMMARY STATISTICS ROWNAMES IN DATA FOR FUTURE SEPARATION

htseq_rows <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
match(htseq_rows, rownames(annot)) -> noncount_rows

# REMOVE HTSEQ-COUNT SUMMARY STATISTICS ROWS FROM THE ANNOTATION TABLE

annot[-c(noncount_rows),] -> annot


# FOR EACH HTSEQ-COUNT FILE, READ COUNTS

sapply(htseqcount_files, function(x){
  
  counts <- read.delim(x, header = F,sep = "\t", stringsAsFactors = F,row.names = 1)[,4]
  counts_only <- counts[-c(noncount_rows)]

  return(counts_only)
}) -> count_matrix

# FOR EACH HTSEQ-COUNT FILE, READ SUMMARY STATISTICS

sapply(htseqcount_files, function(x){
  
  counts <- read.delim(x, header = F,sep = "\t", stringsAsFactors = F,row.names = 1)[,4]
  noncount_data <- counts[noncount_rows]
  
  return(noncount_data)
}) -> noncount_matrix

# ADD ANNOTATION INFORMATION TO THE COUNT MATRIX AND SET SOME COLUMN NAMES

cbind(annot, count_matrix) -> count_matrix_with_annot

colnames(count_matrix_with_annot)[1] <- "Gene_Name"
colnames(count_matrix_with_annot)[2] <- "Gene_Symbol"
colnames(count_matrix_with_annot)[3] <- "Gene_Description"

# ADD ANNOTATION INFORMATION TO SUMMARY STATISTICS

cbind(htseq_rows, noncount_matrix) -> htseq_summary_matrix

# WRITE RAW COUNT FILES, SUMMARY STATISTICS

write.table(count_matrix_with_annot, file=paste("raw.counts.",opt$outfile,sep=""), sep="\t", quote=F, col.names = NA)
write.table(htseq_summary_matrix, file=paste("htseq-seq_summary.",opt$outfile,sep=""), sep="\t", quote=F, col.names = NA)

# LOAD SAMPLE INFORMATION FILE

sample_info <- read.delim(opt$sampleinfofile, header = T,stringsAsFactors = F,row.names=1)
sample_info <- sample_info[colnames(count_matrix),]

# LOADING DATA INTO DESEQ2 COUNT MATRIX OBJECT

rownames(count_matrix) <- rownames(annot)
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



