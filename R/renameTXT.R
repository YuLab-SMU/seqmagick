##' rename txt file (eg Description line of fasta file)
##' according to first token (eg accession number)
##'
##'
##' @title renameTXT
##' @param txt_file txt file
##' @param name_file name file
##' @param sep separator
##' @param split logical, split result or not
##' @return None
##' @export
##' @author Guangchuang Yu
renameTXT <- function(txt_file, name_file, sep="_", split=TRUE) {
    txt <- readLines(txt_file)
    res <- renameTXT_internal(txt, name_file, sep)
    if (split) {
        f1 <- file(gsub("\\.", "_IN\\.", txt_file), "w")
        writeLines(res$name[res$ii], f1)
        close(f1)
        f2 <- file(gsub("\\.", "_OUT\\.", txt_file), "w")
        writeLines(res$name[res$jj], f2)
        close(f2)
    } else {
        outfile <- file(gsub("\\.", "_rename\\.", txt_file), "w")
        writeLines(res$name, outfile)
        close(outfile)
    }
}

renameTXT_internal <- function(txt, name_file, sep="_") {
    sep <- paste0("\\", sep)
    acc_txt <- get_first_token(gsub("^>", "", txt), sep)

    nm <- readLines(name_file)
    nm <- unique(nm)
    acc <- get_first_token(gsub("^>", "", nm), sep)

    ii <- jj <- c()
    for (i in seq_along(acc_txt)) {
        j <- which(acc_txt[i] == acc)
        if (length(j) == 1) {
            ii <- c(ii, i)
            txt[i] <- nm[j]
        } else if (length(j) > 1) {
            ii <- c(ii, i)
            txt[i] <- nm[j[1]]
            warning("multiple match, first record will be used...\n", paste0("\t", nm[j], "\n"))
        } else {
            jj <- c(jj, i)
        }
    }
    res <- list(name=txt, ii=ii, jj=jj)
    return(res)
}

get_first_token <- function(x, sep) {
    trimws(sapply(x, function(y) unlist(strsplit(y, sep))[1]))
}
