##' convert phylip file to fasta file
##'
##' 
##' @title phy2fas
##' @param phyfile phylip file 
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return NULL
##' @export
##' @author Guangchuang Yu
phy2fas <- function(phyfile, outfile="out.fas", type="interleaved") {
    x <- phy_read(phyfile)
    fa_write(x, outfile, type)
}
