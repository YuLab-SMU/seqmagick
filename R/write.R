##' write fasta file
##'
##' 
##' @title fa_write
##' @param x XStringSet object
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return NULL
##' @importFrom Biostrings writeXStringSet
##' @export
##' @author Guangchuang Yu
##' @references \url{http://evolution.genetics.washington.edu/phylip/doc/sequence.html}
fa_write <- function(x, outfile, type="interleaved") {
    seq_write(x, filepath=outfile, type, format="fasta")
}

seq_write <- function(x, filepath, type="interleaved", format) {
    type <- match.arg(type, c("interleaved", "sequential"))
    if (type == "interleaved") {
        .write <- writeXStringSet
    } else {
        .write <- writeXStringSet2
    }
    .write(x, filepath, format = format)
}

writeXStringSet2 <- function(x, filepath, format="fasta") {
    out <- file(filepath, "w")

    desc <- names(x)

    if (format == "fasta") {
        for (i in seq_along(desc)) {
            writeLines(paste('>', desc[i]), out)
            writeLines(toString(x[i]), out)
        }
    }
    
    close(out)
}

