##' replace character for example from '-' to 'N' of fasta sequence that only applied inside sequence
##' any '-' character at start/end of the sequence (aligned seqs may contains '-' at prefix/suffix) will not be replaced
##' 
##'
##' 
##' @title replaceInside
##' @param fasfile fasta file
##' @param from character to be replaced, '-' by default
##' @param to replace 'from' to 'to', 'N' by default
##' @param outfile output file
##' @importFrom Biostrings readDNAStringSet
##' @importFrom Biostrings writeXStringSet
##' @export
##' @return DNAStringSet
##' @author Guangchuang Yu
replaceInside <- function(fasfile, from="-", to="N", outfile=NULL) {
    fas <- readDNAStringSet(fasfile)
    fas_replaced <- replaceInside.DNAStringSet(fas, from, to)
    if (is.null(outfile)) {
        outfile <- sub("\\.(fa.*)$", "_replaced.\\1", fasfile)
    }
    writeXStringSet(fas_replaced, outfile)
    invisible(fas_replaced)
}

##' @importFrom Biostrings DNAStringSet
##' @importFrom magrittr %>%
replaceInside.DNAStringSet <- function(x, from, to) {
    strings <- sapply(1:length(x), function(i) toString(x[i])) %>%
        sapply(., replaceInside.character, from=from, to=to)
    names(strings) <- names(x)
    DNAStringSet(strings)
}


replaceInside.character <- function(x, from="-", to="N") {
    ## pattern matches start/end will be ignore
    ##
    ## > replaceInside.character("----ABC--ABDC---EDFFF---")
    ## [1] "----ABCNNABDCNNNEDFFF---"
    
    y <- gsub(from, to, x)
    
    x <- gsub(paste0("^", to), from, y)

    ## restore start
    while (x != y) {
        y <- x
        x <- gsub(paste0("^(", from, "*)", to), paste0("\\1", from), y)
    }

    ## restore end
    y <- x
    x <- gsub(paste0(to, "$"), from, y)
    while(x != y) {
        y <- x
        x <- gsub(paste0(to, "(", from, "*)$"), paste0(from, "\\1"), y)
    }

    return(x)
}
