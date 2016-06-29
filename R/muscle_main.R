#####################################################################################
# R wrapper function(s) for Edgar's multiple sequence alignment (MUSCLE) software.  #
# see:                                                                              #
# Edgar, R.C. (2004)                                                                #
#	MUSCLE: multiple sequence alignment with high accuracy and high throughput. #
#	Nucleic Acids Res 32, 1792-1797.                                            #
# and:                                                                              #
# http://www.drive5.com/muscle/muscle_userguide3.8.html                             #
#                                                                                   #
# Original author of MUSCLE algorithm: Robert C. Edgar                              #
# Ported into R by Alex T. Kalinka (alex.t.kalinka@gmail.com)                       #
#                                                                                   #
#####################################################################################


###
### from https://raw.githubusercontent.com/Bioconductor-mirror/muscle/master/R/muscle_main.R
###
### modified by Guangchuang Yu to support input stringset as BStringSet
###

##' @importFrom Biostrings readAAStringSet
##' @importFrom Biostrings readDNAStringSet
##' @importFrom Biostrings readRNAStringSet
##' @importFrom Biostrings AAMultipleAlignment
##' @importFrom Biostrings DNAMultipleAlignment
##' @importFrom Biostrings RNAMultipleAlignment
muscle2 <- function(stringset, quiet = FALSE, ...) {
    ## stringset is an object of class 'XStringSet': DNAStringSet, RNAStringSet, or AAStringSet.
    
    args <- list(...)
    arg.v <- NULL
    
    cl <- class(stringset)
    if(length(grep("StringSet", cl))==0){
        stop("input must be an object of class XStringSet: DNAStringSet, RNAStringSet, or AAStringSet\n")
    }
    
    ## Get names and sequence order for re-ordering post-algorithm.
    snames <- names(stringset)
    names(stringset) <- 1:length(stringset)
    
    tempIn <- tempfile(fileext = ".fa")
    writeXStringSet(stringset, filepath = tempIn, format = "fasta")
    arg.v <- append(arg.v,c("in", as.character(tempIn)))
    
    tempOut <- tempfile(fileext = ".afa")
    arg.v <- append(arg.v, c("out", tempOut))
    
    if (quiet){
        arg.v <- append(arg.v, as.character("quiet"))
    }
    if(length(args) > 0){
        for(i in 1:length(args)){
            if(is.logical(args[[i]])){
                ## Is a flag.
                arg.v <- append(arg.v, as.character(names(args)[i]))
            } else{ 
                ## Is an option.
                arg.v <- append(arg.v, c(as.character(names(args)[i]),as.character(args[[i]])))
            }
        }
    }
    nargs <- as.integer(length(arg.v))

    stuff <- .C(.muscleR, PACKAGE="muscle",
                nargs, as.character(arg.v))
     
    ret <- switch(
        class(stringset),
        BStringSet = readBStringSet(tempOut, format = "fasta"),
        AAStringSet = readAAStringSet(tempOut, format = "fasta"),
        DNAStringSet = readDNAStringSet(tempOut, format = "fasta"),
        RNAStringSet = readRNAStringSet(tempOut, format = "fasta")
    )
    
    ## Re-order and re-name output.
    ret <- ret[order(as.integer(names(ret)))]
    names(ret) <- snames

    ret <- switch(
        class(stringset),
        BStringSet = ret,
        AAStringSet = AAMultipleAlignment(ret),
        DNAStringSet = DNAMultipleAlignment(ret),
        RNAStringSet = RNAMultipleAlignment(ret)
    )

    file.remove(tempIn, tempOut)
    
    return(ret)
    
}




