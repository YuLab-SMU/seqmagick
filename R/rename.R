##' rename fasta sequence name
##'
##' 
##' @title fa_rename
##' @param fasfile fasta file
##' @param mapping_file mapping file
##' @param position rename token at specific position
##' @param sep sepator to divide token
##' @param replace logical, whether replace token or not
##' @param prefix logical. If replace is FALSE, whether add to prefix or suffix
##' @param outfile output file
##' @return BStringSet object
##' @importFrom utils read.delim
##' @export
##' @author Guangchuang Yu
fa_rename <- function(fasfile, mapping_file, position, sep, replace, prefix, outfile) {
    x <- fa_read(fasfile)
    mapping <- read.delim(mapping_file, header=FALSE, stringsAsFactors=FALSE)
    y <- bs_rename(x, mapping, position, sep, replace, prefix)
    if (!missing(outfile)) {
        fa_write(y, outfile)
    }
    invisible(y)
}

##' rename sequence
##'
##' 
##' @title bs_rename
##' @param x BStringSet object
##' @param mapping two column data.frame
##' @param position rename token at specific position
##' @param sep sepator to divide token
##' @param replace logical, whether replace token or not
##' @param prefix logical. If replace is FALSE, whether add to prefix or suffix
##' @return BStringSet
##' @export
##' @author Guangchuang Yu
bs_rename <- function(x, mapping, position, sep, replace, prefix) {
    nn <- names(x)
    
    if (missing(sep)) {
        id <- names(x)
    } else {
        id <- sapply(tokenize(nn, sep), function(y) y[1])
    }
    
    if (!missing(position)) {
        if (missing(sep))
            stop("separator should be provided")
        
        tk <- sapply(tokenize(nn, sep), function(y) y[position])
        
        i <- match(id, mapping[,1])
        
        if (replace) {
            to <- mapping[i,2]
        } else if (prefix) {
            to <- paste(mapping[i, 2], tk, sep="")
        } else {
            to <- paste(tk, mapping[i, 2], sep="")
        }
        
        tks <- tokenize(nn, sep)
        tks <- lapply(seq_along(tks), function(i) {
            tks[[i]][position] <- to[i]
            tks[[i]]
        })

        id2 <- sapply(tks, paste0, collapse=sep)
        
    } else {
        if (replace) {
            i <- match(id, mapping[,1])
            id2 <- mapping[i, 2]
        } else if (prefix) {
            id2 <- paste(mapping[i,2], id, sep="")
        } else {
            id2 <- paste(id, mapping[i, 2], sep="")
        
        }
    }
    
    names(x) <- id2
    return(x)
}

tokenize <- function(x, sep=" ") {
    lapply(x, function(y) unlist(strsplit(y, split=sep)))
}
