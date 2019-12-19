##' write fasta file
##'
##' 
##' @title fa_write
##' @param x XStringSet object
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return None
##' @export
##' @author Guangchuang Yu
##' @references \url{http://evolution.genetics.washington.edu/phylip/doc/sequence.html}
##' @examples
##' phy_file <- system.file("extdata/HA.phy", package="seqmagick")
##' x <- phy_read(phy_file)
##' fa_file <- tempfile(fileext = '.fas')
##' fa_write(x, fa_file)
fa_write <- function(x, outfile, type="interleaved") {
    type <- match.arg(type, c("interleaved", "sequential"))
    if (type == "interleaved") {
        fa_write_interleaved(x, outfile)
    } else {
        fa_write_sequential(x, outfile)
    }
}

##' write phylip file
##'
##' 
##' @title phy_write
##' @param x XStringSet object
##' @param outfile output file
##' @param type one of interleaved and sequential 
##' @return None
##' @export
##' @author Guangchuang Yu
##' @examples
##' \dontrun{
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' aln <- bs_aln(x[1:5])
##' phy_file <- tempfile(fileext = '.phy')
##' phy_write(aln, phy_file)
##' }
phy_write <- function(x, outfile, type="sequential") {
    type <- match.arg(type, c("interleaved", "sequential"))

    if (type == "interleaved") {
        phy_write_interleaved(x, outfile)
    } else {
        phy_write_sequential(x, outfile)
    }
}

check_alnseq <- function(x) {
    if (length(unique(width(x))) != 1) {
        stop("input fasfile should contains aligned sequences with equal length...")
    }
}

##' @importFrom Biostrings width
##' @importFrom Biostrings toString
phy_write_sequential <- function(x, outfile) {
    check_alnseq(x)

    w <- width(x[1])
    out <- file(outfile, "w")
    header <- paste(length(x), "\t", w)
    writeLines(header, out)
    nn <- max(nchar(names(x)))
    for (i in 1:length(x)) {
        n <- names(x[i])
        sep.blank <- paste(rep(" ", nn-nchar(n)+4), sep="", collapse="")
        line <- paste(n, sep.blank, toString(x[i]), sep="")
        writeLines(line, out)
    }
    close(out)
}

phy_write_interleaved <- function(x, outfile) {
    check_alnseq(x)
    
    w <- width(x[1])
    out <- file(outfile, "w")
    header <- paste(length(x), "\t", w)
    writeLines(header, out)
    nn <- max(nchar(names(x)))
    for (j in 1:ceiling(w/50)) {
        for (i in seq_along(x)) {
            if (j == 1) {
                n <- names(x[i])
                sep.blank <- paste(rep(" ", nn-nchar(n)+4), sep="", collapse="")
                preceding <- paste(n, sep.blank)
            } else {
                preceding <- paste(rep(" ", nn + 5), sep="", collapse="")
            }
            
            to <- seq(j*5-4, j*5) * 10
            from <- to - 9
            if (to[5] > w) {
                idx <- which(to >= w)[1]
                from <- from[1:idx]
                to <- to[1:idx]
                to[idx] <- w
            }
            s <- paste(substring(toString(x[i]), from, to), collapse=" ")
            line <- paste0(preceding, s)
            writeLines(line, out)
        }
        if (j != ceiling(w/50)) {
            linebreak=""
            writeLines(linebreak, out)
        }
    }
    close(out)
}

##' @importFrom Biostrings writeXStringSet
fa_write_interleaved <- function(x, outfile) {
    writeXStringSet(x, filepath=outfile)
}

fa_write_sequential <- function(x, outfile) {
    out <- file(outfile, "w")
    desc <- names(x)

    for (i in seq_along(desc)) {
        writeLines(paste('>', desc[i]), out)
        writeLines(toString(x[i]), out)
    }
    
    close(out)
}

