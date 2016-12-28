##' convert bam file to aligned fasta file
##'
##'
##' @title bam2DNAStringSet
##' @param bamfile bam file
##' @param refseq refseq, DNAStringSet object
##' @return DNAStringSet object
##' @importFrom GenomicAlignments stackStringsFromBam
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges IRanges
##' @importFrom Rsamtools BamFile
##' @export
##' @author guangchuang yu
bam2DNAStringSet <- function(bamfile, refseq) {
    w <- width(refseq)
    if (length(refseq) > 1) {
        stop("please use bam2DNAStringSet2...")
    }
    pm <- GRanges(names(refseq), IRanges(1, w))

    stackStringsFromBam(BamFile(bamfile), param=pm,
                        Lpadding.letter='-', Rpadding.letter='-', use.names=TRUE)

}

##' convert bam file to aligned fasta file
##'
##'
##' @title bam2DNAStringSet2
##' @param bamfile bam file
##' @param refseq refseq, DNAStringSet object
##' @return DNAStringSet object
##' @export
##' @author guangchuang yu
bam2DNAStringSet2 <- function(bamfile, refseq) {
    lapply(seq_along(refseq), function(i) {
        bam2DNAStringSet(bamfile, refseq[i])
    })
}
