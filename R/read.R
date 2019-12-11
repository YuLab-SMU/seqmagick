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
    a <- strsplit(toupper(string), split="")[[1]]
    freq <- mean(a %in% c('A', 'C', 'G', 'T', 'X', 'N', '-') )
    if (freq > 0.9) {
        return('DNA')
    }

    freq <- mean(a %in% c('A', 'C', 'G', 'U', 'X', 'N', '-') )
    if (freq > 0.9) {
        return('RNA')
    }

    return('AA')
}
