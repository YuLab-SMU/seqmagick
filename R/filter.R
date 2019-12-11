##' fasta filter by searching pattern
##'
##' 
##' @title fa_filter
##' @param fasfile input fasta file
##' @param pattern keyword for filter
##' @param by one of 'description' and 'sequence'
##' @param ignore.case logical
##' @param outfile output file
##' @param type one of 'interleaved' and 'sequential'
##' @return BStringSet object
##' @importFrom Biostrings readBStringSet
##' @export
##' @author Guangchuang Yu
fa_filter <- function(fasfile, pattern, by='description', ignore.case=FALSE,
                      outfile=NULL, type="interleaved") {
    
    fas <- readBStringSet(fasfile)

    res <- bs_filter(fas, pattern, by, ignore.case)
    
    if (!is.null(outfile)) {
       fa_write(res, outfile, type)
    }
    
    invisible(res)
}

##' biological sequence filter by searching pattern
##'
##' 
##' @title bs_filter
##' @param x BStringSet object
##' @param pattern keyword for filter
##' @param by one of 'description' and 'sequence'
##' @param ignore.case logical
##' @return BStringSet object
##' @importFrom Biostrings matchPattern
##' @importFrom magrittr %<>%
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' bs_filter(x, 'ATGAAAGTAAAA', by='sequence')
bs_filter <- function(x, pattern, by='description', ignore.case=FALSE) {
    by <- match.arg(by, c("description", "sequence"))
    
    if (by == 'description') {
        desc <- names(x)
        if (ignore.case) {
            pattern %<>% toupper
            desc %<>% toupper
        }
        idx <- grep(pattern, desc)
    } else {
        if (ignore.case) {
            hit <- sapply(x, function(y) matchPattern(toupper(pattern), toupper(y)))
        } else {
            hit <- sapply(x, function(y) matchPattern(pattern, y))
        }
        idx <- which(sapply(hit, length) > 0)
    }

    if (length(idx) == 0) {
        message("pattern not found...")
        return(NULL)
    }
    
    return(x[idx])
}
