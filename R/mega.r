#' @rdname msa-read
#' @export
#' @importFrom yulab.utils use_perl
mega_read <- function(file) {
    x <- readLines(file)
    start <- grep("^#", x, perl=use_perl())[2] - 1
    header <- x[1:start]
    format <- mega_format(header)

    if (is.null(format)) {
        res <- mega_record_no_format(x[-(1:start)])
        return(res)
    }

    i <- grep("^!Gene=", x, perl=use_perl())

    if (length(i) == 0) {
        records <- x[-(1:start)]

        res <- mega_record(records, format)
        return(res)
    } 
    genename <- sub("!Gene=", "", x[i], perl=use_perl())
    genename <- sub(";", "", genename, perl=use_perl())
    genename <- sub("\\s+.*", "", genename, perl=use_perl())

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
    records <-  sub("^#", ">", records, perl=use_perl())
    records <- gsub("\\s+", "", records, perl=use_perl())

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
    x <- parse_name_seq(records)
    seqname <- sub("^#", "", x$seqname, perl=use_perl())
    seqs <- collapse_seqs(seqname, x$seqs)
    
    build_seq_object(seqs, type)
}

collapse_seqs <- function(seqname, seqs) {
    if (any(duplicated(seqname))) {
        seqs <- vapply(split(seqs, seqname), paste0, collapse = "", FUN.VALUE=character(1))
    } else {
        names(seqs) <- seqname
    }
    return(seqs)
}

parse_name_seq <- function(records) {
    records <- records[records != ""]
    seqname <- sub("^(\\S+)\\s+.*", "\\1", records, perl=use_perl())
    seqs <- sub("^\\S+", "", records, perl=use_perl())
    seqs <- gsub("\\s+", "", seqs, perl=use_perl())
    return(list(seqname = seqname, seqs = seqs))
}


build_seq_object <- function(seq_str, type) {
    if (is.na(type) || type == "auto") {
        type <- guess_sequence_type(seq_str[1])
    }

    type <- toupper(type)
    
    switch(type,
           DNA = DNAStringSet(seq_str),
           RNA = RNAStringSet(seq_str),
           AA  = AAStringSet(seq_str),
           PROTEIN  = AAStringSet(seq_str),
           UNKNOWN = BStringSet(seq_str)
           )    
}

mega_format <- function(header) {
    i <- grep("!Format\\s+DataType=", header, perl=use_perl())
    if (length(i) == 0) {
        # not format definition
        return(NULL)
    } 

    f <- gsub("!Format\\s*|;", "", header[i], perl=use_perl())
    ff <- strsplit(strsplit(f, " ")[[1]], "=") 
    fm <- do.call(rbind, ff)
    res <- setNames(fm[,2], fm[,1])
    return(res)
}
