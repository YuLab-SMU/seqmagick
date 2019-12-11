##' hamming distances of sequences
##'
##'
##' @title bs_hamming
##' @param x BStringSet object
##' @param count_indel whether count indel or not
##' @param ... additional parameters
##' @return hamming distance
##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' ## align first 5 sequences, use `bs_aln(x)` to align all sequences
##' aln <- bs_aln(x[1:5])
##' bs_hamming(aln)
bs_hamming <- function(x, count_indel=FALSE, ...) {
    n <- length(x)
    if (n < 2) {
        stop("at least 2 sequences was expected...")
    }

    if ( n == 2 ) {
        s1 <- toString(x[1])
        s2 <- toString(x[2])
        res <- str_hamming(s1, s2, count_indel, ...)
        return(res)
    }

    res <- matrix(NA, nrow=n, ncol=n)
    colnames(res) <- rownames(res) <- names(x)
    seqs <- lapply(x, function(xx) {
        y <- toString(xx)
        substring(y, 1:nchar(y), 1:nchar(y))
    })
    cnt <- 1
    pb <- txtProgressBar(min=0, max=sum(1:n), style=3)

    for (i in 1:n) {
        for (j in 1:i) {
            setTxtProgressBar(pb, cnt)
            cnt <- cnt + 1
            if (i == j) {
                res[i,j] <- 0
            } else {
                res[i,j] <- suppressMessages(charVector_hamming(seqs[[i]], seqs[[j]], count_indel, ...))
                res[j, i] <- res[i, j]
            }
        }
    }
    close(pb)
    return(res)
}

str_hamming <- function(seq1, seq2, count_indel=FALSE, ...) {
    if (length(seq1) == 1 & length(seq2) == 1) {
        n1 <- nchar(seq1)
        n2 <- nchar(seq2)
        if (n1 != n2) {
            stop("length of two sequences should be consistent...")
        }
        seq1 <- substring(seq1, 1:n1, 1:n1)
        seq2 <- substring(seq2, 1:n2, 1:n2)
    }

    if ( length(seq1) > 1 & length(seq2) > 1) {
        charVector_hamming(seq1, seq2, count_indel, ...)
    } else {
        stop("seq1 and seq2 should be character of equal length...")
    }
}


charVector_hamming <- function(seq1, seq2, count_indel, type = "mapping", threshold=0.3) {
    type <- match.arg(type, c("mapping", "alignment"))
    ## mapping for mapping short read to long sequence
    ## alignment for aligning sequences of similar lengths

    ## for mapping
    ##
    ##
    ## gap will not be counted
    ##
    ## ---------XXXXXXX----------
    ## XXXXXXXXXXXXXXXXXXXXXXXXXX
    ##
    ## gap should be counted
    ##
    ## -------------------XXXXXXX
    ## XXXXXXXXXXXXXXXXXXXXX-----
    ##

    if (type == "mapping") {
        n1 <- seq1 != '-'
        n2 <- seq2 != '-'
        nn <- min(sum(n1), sum(n2))
        if (sum(n1 & n2)/nn < threshold) {
            count_indel <- TRUE
        }
    }

    ii <- which(seq1 != seq2)

    if (count_indel == FALSE) {
        message("--> indel will not count...")
        ii <- ii[seq1[ii] != '-' & seq2[ii] != '-']
    }

    return(length(ii))
}

