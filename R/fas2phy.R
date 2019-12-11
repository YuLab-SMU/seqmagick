##' convert fasta (aligned sequences) to phylip format
##'
##' 
##' @title fas2phy 
##' @param fasfile aligned sequences in fasta format
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return None
##' @export
##' @author Guangchuang Yu
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' phy_file <- tempfile(fileext = ".phy")
##' fas2phy(fa_file, phy_file)
fas2phy <- function(fasfile, outfile="out.phy", type="sequential") {
    x <- fa_read(fasfile)
    phy_write(x, outfile, type)
}
