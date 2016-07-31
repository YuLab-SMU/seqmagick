##' read fasta file
##'
##' 
##' @title fa_read
##' @param file fasta file
##' @return BStringSet object
##' @importFrom Biostrings readBStringSet
##' @export
##' @author Guangchuang Yu
fa_read <- function(file) {
    readBStringSet(file)
}


##' read aligned sequences in phylip format
##'
##' 
##' @title phy_read
##' @param file phylip file
##' @importFrom Biostrings BStringSet
##' @return BStringSet object
##' @export
##' @author ygc
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
        seq_str <- substring(phy, ii+1, info$widt)
        
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
    BStringSet(seq_str)
}

