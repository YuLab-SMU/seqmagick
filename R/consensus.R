##' consensus of aligned sequences
##'
##' 
##' @title consensus
##' @rdname consensus
##' @param x XStringSet or XMultipleAlignment object
##' @param type currently, only DNA supported
##' @return consensus sequence string
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' ## align first 5 sequences, use `bs_aln(x)` to align all sequences
##' aln <- bs_aln(x[1:5])
##' ## or bs_consensus(aln)
##' consensus(aln)
consensus <- function(x, type="DNA") {
    BStringSet(x) %>% bs_consensus
}

##' consensus of aligned sequences
##'
##' 
##' @title bs_consensus
##' @rdname consensus
##' @param x BStringSet object
##' @param type currently, only DNA supported
##' @param r if any NT > r, it will be selected as representative base
##' @return consensus sequence string
##' @importFrom Biostrings width
##' @importFrom Biostrings consensusMatrix
##' @export
##' @author Guangchuang Yu
bs_consensus <- function(x, type="DNA", r=1) {
    if (length(unique(width(x))) != 1) {
        stop("input sequences were not aligned...")
    }
    
    if (type != "DNA")
        stop("currently only DNA is supported...")
    
    y <- consensusMatrix(x)
    cns <- rep(NA, ncol(y))

    if (type == "DNA") {
        alphabet <- c("A", "C", "G", "T")
        ambiguity_code <- DNA_ambiguity_code
    }
    
    idx <- match(alphabet, rownames(y))
    nn <- apply(y[idx,], 2, function(n) sum(n != 0))
    i <- which(nn == 1)
    j <- apply(y[idx, i], 2, function(n) which(n!=0))
    cns[i] <- alphabet[j]
    
    cns[nn == 0] <- '-'
    
    if (sum(nn > 1)> 0) {
        ambiguity_alphabet <- apply(y[1:4, nn > 1, drop=FALSE], 2, function(n) {
            rr <- n/sum(n)
            if (any(rr > r)) {
                ii <- which(rr > r)
            } else {
                ii <- which(rr > 0)
            }
            ambiguity_code(alphabet[ii])
        })
        cns[nn > 1] <- ambiguity_alphabet
    }
    paste0(cns, collapse='')
}

##' @importFrom magrittr %>%
DNA_ambiguity_code <- function(NT) {
    ## M          A or C
    ## R          A or G
    ## W          A or T
    ## S          C or G
    ## Y          C or T
    ## K          G or T
    ## V        A or C or G
    ## H        A or C or T
    ## D        A or G or T
    ## B        C or G or T
    ## X      G or A or T or C
    ## N      G or A or T or C

    NT <- NT %>% unique %>% toupper %>% sort
    
    if (length(NT) < 2) {
        return(NT)
    }

    if (length(NT) > 4) {
        stop("length of 'NT' should be ranged between 2 and 4...")
    }

    if (length(NT) == 2) {
        if (all(NT == c("A", "C")))
            return("M")
        if (all(NT == c("A", "G")))
            return("R")
        if (all(NT == c("A", "T")))
            return("W")
        if (all(NT == c("C", "G")))
            return("S")
        if (all(NT == c("C", "T")))
            return("Y")
        if (all(NT == c("G", "T")))
            return("K")
    } else if (length(NT) == 3) {
        if (all(NT == c("A", "C", "G")))
            return("V")
        if (all(NT == c("A", "C", "T")))
            return("H")
        if (all(NT == c("A", "G", "T")))
            return("D")
        if (all(NT == c("C", "G", "T")))
            return("B")
    } else {
        if (all(NT == c("A", "C", "G", "T")))
            return("N")
        
    }
    
    stop("wrong alphabet...")
}


