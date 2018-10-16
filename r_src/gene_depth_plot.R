# plot gene region depth

args<-commandArgs(T)
if (length(args) != 6){
  print('Rscript this.R <infile.csv> <contig> <start> <end> <bin> <outprefix>')
  q()
}

file = args[1]
contig = as.character(args[2])
starter =as.numeric(args[3])
ender=as.numeric(args[4])
bin=as.numeric(args[5])
output=args[6]


library(Gviz)
data(geneModels) 

gtr <- GenomeAxisTrack()
itr <- IdeogramTrack(genome="hg19", chromosome=contig)
grtr<-BiomartGeneRegionTrack(genome="hg19",
                             chromosome=contig,name="Ensembl",
                             from=starter,end=ender)

gr<-GRanges(seqnames = contig, 
           ranges = IRanges(start = seq(starter,ender-bin,by=bin), 
                            width = bin))
mdt = read.csv(file)
mdt=mdt[-1:-3]
values(gr) = mdt


dTrack <- DataTrack(gr, name="depth",type = c('a','heatmap'))

pdf(paste(output,'geneTrackDepth.pdf',sep=''))
plotTracks(list(itr,gtr,dTrack,grtr),from=starter-4000,to=ender+4000,showId=TRUE,showBandId=TRUE)
dev.off()


