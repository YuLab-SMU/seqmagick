
getPhyInfo <- function(phy) {
    line1 <- readLines(phy, n=1)
    res <- strsplit(line1, split="\\s")[[1]]
    res <- res[res != ""]

    return(list(num=as.numeric(res[1]), width=as.numeric(res[2])))
}
