##' download genbank or fasta file by accession number
##'
##' 
##' @title download_genbank
##' @param acc accession number(s)
##' @param db supported db, currently 'nuccore'
##' @param format one of 'genbank' or 'fasta'
##' @param outfile output file, by default, acc.gb or acc.fa
##' @return output file vector
##' @importFrom magrittr %>%
##' @export
##' @author Guangchuang Yu
download_genbank <- function(acc, db="nuccore", format = "genbank", outfile=NULL) {
    if (!is.null(outfile) && (length(outfile) != length(acc))) {
        stop("'outfile' & 'acc' should be equal length...")
    }
    acc <- gsub("\\.\\d+", "", acc)

    ofs <- character(length(acc))
    for (i in seq_along(acc)) {
        if (is.null(outfile)) {
            of <- NULL
        } else {
            of <- outfile[i]
        }
        ofs[i] <- download_genbank_item(acc[i], db, format, of)
    }
    invisible(ofs)
}

##' @importFrom utils download.file
download_genbank_item <- function(acc, db="nuccore", format = "genbank", outfile=NULL) {
    format <- match.arg(format, c('genbank', 'fasta'))
    if (db != "nuccore") {
        stop("currently, only nuccore is supported...")
    }
    url <- paste0("www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=",
                  db, '&dopt=', format, '&sort=&val=', acc)
    
    if (is.null(outfile)) {
        suffix <- '.gb'
        if (format == 'fasta')
            suffix <- '.fa'
        outfile <- paste0(acc, suffix)
    }
    
    download.file(url=url, destfile = outfile, method="curl")
    invisible(outfile)
}

