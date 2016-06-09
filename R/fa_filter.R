##' fasta filter by search pattern
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
##' @importFrom Biostrings matchPattern
##' @importFrom Biostrings writeXStringSet
##' @importFrom magrittr %<>%
##' @export
##' @author Guangchuang Yu
fa_filter <- function(fasfile, pattern, by='description', ignore.case=FALSE,
                      outfile=NULL, type="interleaved") {
    
    by <- match.arg(by, c("description", "sequence"))
    
    fas <- readBStringSet(fasfile)

    if (by == 'description') {
        desc <- names(fas)
        if (ignore.case) {
            pattern %<>% toupper
            desc %<>% toupper
        }

        idx <- grep(pattern, desc)
    } else {
        if (ignore.case) {
            hit <- sapply(fas, function(x) matchPattern(toupper(pattern), toupper(x)))
        } else {
            hit <- sapply(fas, function(x) matchPattern(pattern, x))
        }
        idx <- which(sapply(hit, length) > 0)
    }

    if (length(idx) == 0) {
        message("pattern not found...")
        return(NULL)
    }
    
    res <- fas[idx]
    if (!is.null(outfile)) {
       fa_write(res, filepath=outfile, type)
    }
    invisible(res)
}
