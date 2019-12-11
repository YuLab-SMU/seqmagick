##' convert phylip file to fasta file
##'
##' 
##' @title phy2fas
##' @param phyfile phylip file 
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return None
##' @export
##' @author Guangchuang Yu
##' @examples
##' phy_file <- system.file("extdata/HA.phy", package="seqmagick")
##' fa_file <- tempfile(fileext = '.fas')
##' phy2fas(phy_file, fa_file)
phy2fas <- function(phyfile, outfile="out.fas", type="interleaved") {
    x <- phy_read(phyfile)
    fa_write(x, outfile, type)
}
