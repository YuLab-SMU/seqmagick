##' read sequence alignment file from MEGA
##'
##'
##' @title mega_read
##' @param file mega file
##' @return BStringSet object
##' @export
##' @author Guangchuang Yu
##' @examples
##' mega_file <- system.file("extdata/mega/Crab_rRNA.meg", package="seqmagick")
##' mega_read(mega_file)
mega_read <- function(file) {
    x <- readLines(file)
    start <- grep("^#", x)[2] - 1
    header <- x[1:start]
    format <- mega_format(header)

    if (is.null(format)) {
        res <- mega_record_no_format(x[-(1:start)])
        return(res)
    }

    i <- grep("^!Gene=", x)

    if (length(i) == 0) {
        records <- x[-(1:start)]

        res <- mega_record(records, format)
        return(res)
    } 
    genename <- sub("!Gene=", "", x[i])
    genename <- sub(";", "", genename)
    genename <- sub("\\s+.*", "", genename)

    begin <- i + 1
    end <- c(i[-1] -1, length(x))

    res <- lapply(seq_along(genename), function(j) {
        records <- x[begin[j]:end[j]]
        mega_record_no_format(records, format["DataType"])
    })
    names(res) <- genename
    return(res)
}

mega_record <- function(records, format) {
    records <-  sub("^#", ">", records)
    records <- gsub("\\s+", "", records)

    ic <- format['identical']
    type <- format["DataType"]
    if (is.na(type)) {
        type <- guess_sequence_type(records[2])
    }

    f <- tempfile(fileext = ".fa")
    cat(file = f, records, sep="\n")

    if (is.na(ic)) {
        res <- fa_read(f, type)
        return(res)
    }

    bs <- readBStringSet(f)

    bsm <- as.matrix(bs)
    for (i in 2:nrow(bsm)) {
        j = which(bsm[i, ] == ic) 
        bsm[i, j] = bsm[1, j]
    }
    seq_str <- apply(bsm, 1, paste0, collapse='') 

    build_seq_object(seq_str, type)
}

mega_record_no_format <- function(records, type = NA) {
    records <- records[records != ""]
    seqname <- sub("^#(\\S+)\\s+.*", "\\1", records)
    seqs <- sub("^#\\S+", "", records)
    seqs <- gsub("\\s+", "", seqs)
    if (any(duplicated(seqname))) {
        seqs <- vapply(split(seqs, seqname), paste0, collapse = "", FUN.VALUE=character(1))
    } else {
        names(seqs) <- seqname
    }

    if (is.na(type)) {
        type <- guess_sequence_type(seqs[1])
    }
    
    build_seq_object(seqs, type)
}

build_seq_object <- function(seq_str, type) {
    switch(type,
           DNA = DNAStringSet(seq_str),
           RNA = RNAStringSet(seq_str),
           AA  = AAStringSet(seq_str),
           Protein  = AAStringSet(seq_str),
           UNKNOWN = BStringSet(seq_str)
           )    
}

mega_format <- function(header) {
    i <- grep("!Format\\s+DataType=", header)
    if (length(i) == 0) {
        # not format definition
        return(NULL)
    } 

    f <- gsub("!Format\\s*|;", "", header[i])
    ff <- strsplit(strsplit(f, " ")[[1]], "=") 
    fm <- do.call(rbind, ff)
    res <- setNames(fm[,2], fm[,1])
    return(res)
}

