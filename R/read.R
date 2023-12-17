##' read fasta file
##'
##'
##' @title fa_read
##' @param file fasta file
##' @param type one of 'DNA', 'RNA', 'AA', 'unknown' or 'auto'
##' @return BStringSet object
##' @importFrom Biostrings readBStringSet
##' @importFrom Biostrings readDNAStringSet
##' @importFrom Biostrings readRNAStringSet
##' @importFrom Biostrings readAAStringSet
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' fa_read(fa_file)
fa_read <- function(file, type = "auto") {
    type <- match.arg(type, c("DNA", "RNA", "AA", "unknown", "auto" ))
    if (type == "auto") {
        type <- guess_sequence_type(file)
    }
    switch(type,
           DNA = readDNAStringSet(file),
           RNA = readRNAStringSet(file),
           AA  = readAAStringSet(file),
           UNKNOWN = readBStringSet(file)
           )
}


##' read aligned sequences in phylip format
##'
##'
##' @title phy_read
##' @param file phylip file
##' @importFrom Biostrings BStringSet
##' @return BStringSet object
##' @export
##' @author Guangchuang Yu
##' @examples
##' phy_file <- system.file("extdata/HA.phy", package="seqmagick")
##' phy_read(phy_file)
phy_read <- function(file) {
    info <- getPhyInfo(file)
    phy <- readLines(file)

    type <- "interleaved"
    if (length(phy) == (info$num + 1)) {
        type <- "sequential"
    }

    if (type == "sequential") {
        phy <- phy[-1]

        ii <- nchar(phy) - info$width
        nm <- gsub("\\s+$", "", substring(phy, 1,ii-1))
        seq_str <- substring(phy, ii+1, nchar(phy))

        ## may use \t and lines are not in equal length
        ## nm <- gsub("\\s+[^\\s]+$", "", phy)
        ## seq_str <- gsub("^[^\\s]+\\s+", "", phy)

    } else {
        phy <- phy[nchar(phy) != 0]
        phy <- phy[-1]
        ii <- regexpr("\\S", phy[info$num+2])
        nn <- nchar(phy[info$num+2])
        nm <- gsub("\\s+$", "", substring(phy[1:info$num], 1,ii-1))
        ss <- substring(phy, ii, nn)

        seq_str <- sapply(1:info$num, function(i) {
            j <- seq(i, length(ss), by=info$num)
            paste0(ss[j], sep="", collapse="") %>% gsub(" ", "", .)
        })
    }

    names(seq_str) <- nm
    ## BStringSet(seq_str)
    switch(guess_sequence_type(seq_str[1]),
           DNA = DNAStringSet(seq_str),
           RNA = RNAStringSet(seq_str),
           AA = AAStringSet(seq_str))
}

guess_sequence_type <- function(string) {
    if (file.exists(string)){
        seqstr <- readLines(string, n=3)
        if (length(seqstr)==2 || grepl("^>", seqstr[[3]])){
            # >seq1
            # AGCGTACGTGACGTAGCGTAGC
            # >seq2
            a <- seqstr[2]
        }else{
            # >seq1
            # ---------
            # AGCG----C
            # ---------
            # >seq2
            seqstr <- readLines(string, n=20)
            seqind <- grep("^>", seqstr)
            if (length(seqind)==1){
                a <- paste0(seqstr[-1], collapse="")
            }else{
                inds <- seqind[1] + 1
                inde <- seqind[2] - 1
                a <- paste0(seqstr[inds:inde], collapse="")
            }
        }
    }else{
        a <- string
    }
#    a <- strsplit(toupper(a), split="")[[1]]
    guess(a)
#    freq1 <- mean(a %in% c('A', 'C', 'G', 'T', 'X', 'N', '-') )
#    freq2 <- mean(a %in% c('A', 'C', 'G', 'T', 'N'))
#    if (freq1 > 0.9 && freq2 > 0) {
#        return('DNA')
#    }
#
#    freq3 <- mean(a %in% c('A', 'C', 'G', 'U', 'X', 'N', '-'))
#    freq4 <- mean(a %in% c('T'))
#    freq5 <- mean(a %in% c('U'))
#    if (freq3 > 0.9 && freq4==0 && freq5 > 0) {
#        return('RNA')
#    }
#
#    return('AA')
}

guess <- function(seq_str) {
    seq_str <- strsplit(toupper(seq_str), split="")[[1]]
    freq1 <- mean(seq_str %in% c('A', 'C', 'G', 'T', 'X', 'N', '-') )
    freq2 <- mean(seq_str %in% c('A', 'C', 'G', 'T', 'N'))
    if (freq1 > 0.9 && freq2 > 0) {
        return('DNA')
    }
    freq3 <- mean(seq_str %in% c('A', 'C', 'G', 'U', 'X', 'N', '-'))
    freq4 <- mean(seq_str %in% c('T'))
    freq5 <- mean(seq_str %in% c('U'))
    if (freq3 > 0.9 && freq4==0 && freq5 > 0) {
        return('RNA')
    }
    return('AA')
}

##' extract accession number and sequence from genbank file
##'
##' 
##' @title gb_read
##' @param file input genbank file
##' @return sequence object
##' @export
##' @author Guangchuang Yu
gb_read <- function(file) {
    res <- vapply(file, gb_read_item, FUN.VALUE = character(1))
    
    f <- tempfile(fileext = ".fasta")
    cat(res, file = f, sep="")
    fa_read(f)
}

gb_read_item <- function(file) {
    x <- readLines(file)
    acc <- sub("\\w+\\s+(\\w+)$", "\\1", x[grep("^ACCESSION", x)])
    src <- sub("SOURCE\\s+",  "", x[grep("^SOURCE",x)])
    header <- paste0(src, "(", acc, ")")
    i <- grep("ORIGIN", x)
    ss <- x[(i+1):length(x)]
    ss <- ss[1:(grep("//", ss) -1)]
    ss <- gsub("\\s+\\d+", "", ss)
    ss <- gsub("\\s+", "", ss)
    ss <- paste0(ss, collapse = "")

    # f <- tempfile(fileext = ".fasta")
    # cat(">",  file = f, append = TRUE)
    # cat(header,  file = f, append = TRUE, sep = "\n")
    # cat(ss,  file = f, append = TRUE, sep = "\n")
    # fa_read(f)

    sprintf(">%s\n%s\n", header, ss)
}

filt_str <- function(str, to_filt = "") {
    str[str != to_filt]
}

##' read mega file
##'
##'
##' @title mega_read
##' @param file mega file
##' @param type one of 'DNA', 'RNA', 'AA', 'unknown' or 'auto'
##' @return BStringSet object
##' @importFrom Biostrings BStringSet
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings RNAStringSet
##' @importFrom Biostrings AAStringSet
##' @export
##' @author Guangchuang Yu; Zijing Xie
##' @examples
##' mega_file <- system.file("extdata/rRNA.meg", package="seqmagick")
##' mega_read(mega_file)
mega_read <- function(file, type = "auto") {
    mega <- readLines(con = file) |> filt_str()
    mega <- gsub(pattern = "\\s+", replacement = "", mega) |>
        paste(collapse = "\n")
    mega <- strsplit(mega, split = "#")[[1]] |> filt_str()
    mega <- mega[-1]
    raw_result <- strsplit(mega, split = "\n") |> sapply(
        function(x) {
            seq_str <- paste(x[2:length(x)], collapse = "")
            setNames(seq_str, x[1])
        }
    )
    seq_mat <- strsplit(raw_result, split = "") |>
        do.call(cbind, args = _)
    for (i in seq_len(nrow(seq_mat))) {
        seq_mat[i, ] <- gsub(pattern = "\\.",
            replacement = seq_mat[i, 1], seq_mat[i, ])
    }
    seq_str <- apply(seq_mat, 2, function(x) {
        paste(x, collapse = "")
    })
    type <- match.arg(type, c("DNA", "RNA", "AA", "unknown", "auto"))
    if (type == "auto") {
        type <- guess(seq_str[1])
    }
    switch(type,
           DNA = DNAStringSet(seq_str),
           RNA = RNAStringSet(seq_str),
           AA  = AAStringSet(seq_str),
           UNKNOWN = BStringSet(seq_str)
           )
}
