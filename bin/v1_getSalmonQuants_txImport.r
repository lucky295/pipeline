#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
      stop("This script requires at least 2 arguments\nArg1=path to directory that has sample info which then have quant.sf file in
their subfolders\nArg2=Folder path to a directory that has a 2 column file that ends in tx2gene.tsv - column1 is transcript, column2 is gene\n",
call.=FALSE)
}

check_dir_exists<- function(my_path) {
      if (!dir.exists(my_path)){
          cat("This directory does not exist: ", my_path, "\n")
          stop("Exitting - Directory Not Found", call.=FALSE)
      }
}


get_quants <- function(path1, ...) {
    #cat("I am inside get_quants function", path1, "\n")
    additionalPath = list(...)

    suppressMessages(library(tximport))
    suppressMessages(library(readr))

    salmon_filepaths=file.path(path=path1,list.files(path1,recursive=TRUE, pattern="quant.sf"))
    #str(salmon_filepaths)
    samples = data.frame(samples = gsub(".*?quant/salmon_(.*?)/quant.sf", "\\1", salmon_filepaths) )
    row.names(samples)=samples[,1]
    names(salmon_filepaths)=samples$samples
    #print(samples)

    # IF no tx2Gene available, we will only get isoform level counts
    salmon_tx_data = tximport(salmon_filepaths, type="salmon", txOut = TRUE)


    ## Get transcript count summarization
    write.csv(as.data.frame(salmon_tx_data$counts), file = "tx_NumReads.csv")
    ## Get TPM
    write.csv(as.data.frame(salmon_tx_data$abundance), file  =  "tx_TPM_Abundance.csv")
    
    
    if(length(additionalPath > 0)) {
        #cat("tx2gene path available ",  "\n" )
        #print(additionalPath[[1]])
        tx2genePath = additionalPath[[1]]
        tx2geneFile=file.path(path=tx2genePath,list.files(tx2genePath,recursive=FALSE, pattern=".*tx2gene.tsv$"))
        num_quant_tx2geneFile = length(tx2geneFile)

        if (num_quant_tx2geneFile == 1)  {
            #cat("tx2geneFile exists - ", tx2geneFile, sep = "\n")

            my_tx2gene=read.csv(tx2geneFile,sep = "\t",stringsAsFactors = F, header=F)
            salmon_tx2gene_data = tximport(salmon_filepaths, type="salmon", txOut = FALSE, tx2gene=my_tx2gene)
            
            ## Get Gene count summarization
            write.csv(as.data.frame(salmon_tx2gene_data$counts), file = "tx2gene_NumReads.csv")
            ## Get TPM
            write.csv(as.data.frame(salmon_tx2gene_data$abundance),  file  =  "tx2gene_TPM_Abundance.csv")
        }else {
            cat("No filename ending with tx2gene.tsv  found at directory: ", tx2genePath, "\n")
            stop("Exitting", call.=FALSE)
        } 
    }
}

if(length(args) == 1){
    quant_dirpath =   gsub("/$", '', args[1])
    check_dir_exists(quant_dirpath)

    #cat(quant_dirpath, "\n")
    get_quants(quant_dirpath)
}else if(length(args) == 2){

    quant_dirpath =   gsub("/$", '', args[1])
    check_dir_exists(quant_dirpath)

    tx2genePath =  gsub("/$", '', args[2])
    check_dir_exists(tx2genePath)
    
    #cat(quant_dirpath, "\n")
    #cat(tx2genePath, "\n")
    suppressWarnings(get_quants(quant_dirpath, tx2genePath))
    
}else{
    stop("This script requires at least 2 arguments\nArg1=path to directory that has sample info which then have quant.sf file in
      their subfolders\nArg2=Folder path to a directory that has a 2 column file that ends in tx2gene.tsv - column1 is transcript, column2 is gene\n",
call.=FALSE)
}



