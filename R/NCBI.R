##' download genbank or fasta file by accession number
##'
##'
##' @title download_genbank
##' @param acc accession number(s)
##' @param db supported db, currently 'nuccore'
##' @param format one of 'genbank' or 'fasta'
##' @param outfile output file, by default, acc.gb or acc.fa
##' @param ... additional parameters for download.file
##' @return output file vector
##' @importFrom magrittr %>%
##' @export
##' @author Guangchuang Yu
##' @examples
##' \dontrun{
##' tmpgb <- tempfile(fileext = '.gb')
##' tmpfa <- tempfile(fileext = '.fa')
##' download_genbank(acc='AB115403', format='genbank', outfile=tmpgb)
##' download_genbank(acc='AB115403', format='fasta', outfile=tmpfa)
##' }
download_genbank <- function(acc, db="nuccore", format = "genbank", outfile=NULL, ...) {
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
        ofs[i] <- download_genbank_item(acc[i], db, format, of, ...)
    }
    invisible(ofs)
}

##' @importFrom utils download.file
download_genbank_item <- function(acc, db="nuccore", format = "genbank", outfile=NULL, ...) {
    format <- match.arg(format, c('genbank', 'fasta'))
    if (db != "nuccore") {
        stop("currently, only nuccore is supported...")
    }
    url <- paste0("https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=",
                  db, '&dopt=', format, '&sort=&val=', acc)

    if (is.null(outfile)) {
        suffix <- '.gb'
        if (format == 'fasta')
            suffix <- '.fa'
        outfile <- paste0(acc, suffix)
    }
    
    downloader::download(url=url, destfile = outfile,  ...)
    invisible(outfile)
}

##' download fasta sequence by accession number and read as sequence object
##'
##'
##' @title ncbi_fa_read
##' @param acc accession number(s)
##' @param db supported db, currently 'nuccore'
##' @param type one of 'DNA', 'RNA', 'AA', 'unknown' or 'auto'
##' @param ... additional parameters passed to downloader::download
##' @return BStringSet object
##' @export
##' @author Guangchuang Yu
##' @examples
##' \dontrun{
##' x <- ncbi_fa_read(c('AB115403', 'AB115403.1'))
##' x
##' }
ncbi_fa_read <- function(acc, db="nuccore", type = "auto", ...) {
    if (length(acc) == 0) {
        stop("'acc' is empty...")
    }
    if (anyNA(acc)) {
        stop("'acc' contains NA...")
    }

    acc <- as.character(acc)
    acc <- gsub("\\.\\d+", "", acc)

    tmpdir <- tempfile("seqmagick-ncbi-")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive=TRUE, force=TRUE), add=TRUE)

    destfile <- file.path(tmpdir, paste0(make.unique(acc), ".fa"))
    download_genbank(acc, db=db, format="fasta", outfile=destfile, ...)

    f <- tempfile(fileext = ".fasta")
    for (i in seq_along(destfile)) {
        x <- readLines(destfile[i], warn=FALSE)
        if (length(x) == 0) {
            next
        }
        cat(x, file = f, append = TRUE, sep="\n")
        cat("\n", file = f, append = TRUE)
    }

    fa_read(f, type=type)
}

