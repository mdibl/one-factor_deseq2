# simple script for running a single-factor, 2-level experiment using DESeq2 in R from a batch script
# including the use of a json-based configuration file
#
# Load the necessary libraries
suppressMessages(library("rjson"))
suppressMessages(library("DESeq2"))

#get the arguments and create parameters from the json file (expected as the only command-line argment)
args <- commandArgs(trailingOnly = TRUE)
jObj <- fromJSON(file = args[1])
#summarize and print out the complete json object
summary(jObj)
jObj

#read the count data, summarize it, print dimensions, and then round it off 
datafile <- paste(jObj$count_file_directory,jObj$count_file_name,sep="/")
datafile
countData <- as.matrix(read.table(datafile, header = TRUE, sep = jObj$separator, dec = '.', row.names = jObj$gene_col_title,check.names=FALSE))
summary(countData)
dim(countData)
countData <- round(countData)

### change the filtering to the raw data
filter <- apply(countData,1,function(x) mean(x)>jObj$row_sum_th)
table(filter)
countData <- countData[filter,]

#read the design file, print out a summary and the dimensions
designfile <- paste(jObj$design_file_directory,jObj$design_file_name,sep="/")
designfile
colData <- read.table(designfile, header = TRUE, sep = jObj$separator)
summary(colData)
dim(colData)

#create a DESeq data set
dataset <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(paste("",jObj$model_formula, sep = " ~ ")))
#remove rows below the threshold for the sum of all counts across all samples
#dataset <- dataset[rowSums(counts(dataset)) > jObj$row_sum_th, ]
#print a summary
summary(dataset)
dim(dataset)

#Run DESeq
dataset <- DESeq(dataset)
# generate a results object for the contrast listed in the json file
#### Note, this will need to be different for a continuous variable than for a categorical  ####
res <- results(dataset,c(jObj$contrast))
summary(res)

#Dump to an output text (csv) file
write.table(as.data.frame(res),file=paste(jObj$output_directory,paste(jObj$output_prefix,"DESeq2_out.tsv", sep="."),sep="/"),quote=FALSE,sep='\t',col.names=NA,row.names=TRUE)

#print out the normalization factors
sizeFactors(dataset)

#get the dispersions, print a header, and then write to a csv
disp <- dispersions(dataset)
head(disp)
write.table(as.data.frame(disp), file = paste(jObj$output_directory,paste(jObj$output_prefix,"disp.tsv", sep="."),sep="/"),quote=FALSE,sep='\t',col.names=NA,row.names=TRUE)
#perfom, capture, and print out the rlog transform (regularized log2)
rld <- rlog(dataset)
rldAssay <- assay(rld)
write.table(as.data.frame(rldAssay), file = paste(jObj$output_directory,paste(jObj$output_prefix,"rlog.tsv", sep="."),sep="/"),quote=FALSE,sep='\t',col.names=NA,row.names=TRUE)

#generate a file name for pdf pltos of various outputs
pdfname <- paste(jObj$output_directory,paste(jObj$output_prefix,"DESeq2_plots.pdf", sep="."),sep="/")

pdf(pdfname)
plotMA(res)
plotPCA(rld,intgroup=c(jObj$contrast)[1])
plotDispEsts(dataset)
dev.off()

#print out the system/R configuration
sessionInfo()
