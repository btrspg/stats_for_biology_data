args<-commandArgs(T)

if(length(args) != 4){
    print "Rscript coverage_plot.r <coverage.tsv> <names> <chr> <output.pdf>"
    q()
}


library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenVisR")
genomeObject <- BSgenome.Hsapiens.UCSC.hg19
TxDbObject<- TxDb.Hsapiens.UCSC.hg19.knownGene
covData <- read.delim(args[1])

samples = unlist(strsplit(args[2],','))

colnames(covData) <- c("chromosome", "start", "end", samples)

a <- function(x, y){
    col_names <- c("chromosome", "end", x)
    y <- y[,col_names]
    colnames(y) <- c("chromosome", "end", "cov")
    return(y)
}

covData <- lapply(samples, a, covData)

names(covData) <- samples

chromosome <- as.character(unique(covData[[1]]$chromosome))
start <- as.numeric(min(covData[[1]]$end))
end <- as.numeric(max(covData[[1]]$end))

grObject <- GRanges(seqnames=args[3], ranges=IRanges(start=start-500, end=end+500))
pdf(args[4])
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="bar",gene_name="LPL",label_txtSize=7)
dev.off()

