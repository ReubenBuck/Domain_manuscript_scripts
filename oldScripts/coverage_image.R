### correaltion anlysis 




rm(list = ls())


library(GenomicRanges)
library(rtracklayer)
library(zoo)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")
load(file = "~/Desktop/Domain_manuscript/R_objects/chromStateCombined")

### and lets get in the poly morphic ones

polyL1 <- read.csv(file = "~/Desktop/Domain_manuscript/Data/polymorphicL1/polyL1sSIMPLE.csv")
polyL1exact <- polyL1[!is.na(polyL1$Exact.5..Junction),]
polyL1exact <- polyL1exact[,complete.cases(t(polyL1exact))]
polyL1coor <- data.frame(genoName = paste("chr", polyL1exact$Chromosome, sep =""), genoStart = polyL1exact$Exact.5..Junction, genoEnd = polyL1exact$Exact.5..Junction + 6000)

# so lets read in the whole genome and find a way to build the neighbormatrix

# ref gene
web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)

web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chrom_info <- read.delim(textConnection(txt), header = FALSE)

# pull intergenic regions
refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))

refgene_gap.gr <- filterIntergenic(refgene)

intronKeep.gr <- filterIntron(refgene)


# seq gaps
ses <- browserSession("UCSC")
genome(ses) <- genome
data.types <- c("gaps")
track.name <- c("gap")
table.name <- c("gap")
for(i in 1:length(data.types)){
  dat <- getTable(
    ucscTableQuery(
      ses, 
      track = track.name[i],
      table = table.name[i]
    )
  )
  assign(data.types[i], dat)        
}
gaps.gr <- GRanges(seqnames = Rle(gaps$chrom),
                   ranges = IRanges(start = gaps$chromStart, end = gaps$chromEnd)
)		 


# removing gaps intergenic
int <- intersect(refgene_gap.gr, gaps.gr)
OL <- as.matrix(findOverlaps(refgene_gap.gr, int))
agg <- aggregate(x=width(int)[OL[,2]], by=list(OL[,1]), FUN = sum)
bins_gene_gap <- data.frame(chr = seqnames(refgene_gap.gr), start = start(refgene_gap.gr), end = end(refgene_gap.gr), Known = width(refgene_gap.gr))
bins_gene_gap$Known[agg[,1]] = bins_gene_gap$Known[agg[,1]] - agg[,2]
bins_gene_gap <- (bins_gene_gap[bins_gene_gap$Known > bins_gene_gap$end - bins_gene_gap$start -101,])
bins_gene_gap <- (bins_gene_gap[(bins_gene_gap$end - bins_gene_gap$start) / bins_gene_gap$Known > .95,])

# remove gaps from intron
gap.int <- intersect(intronKeep.gr, gaps.gr)
gapOl <- as.matrix(findOverlaps(intronKeep.gr, gap.int))
bins_intron = as.data.frame(intronKeep.gr)[,1:4]
colnames(bins_intron)[c(1,4)] <- c("chr", "Known")
bins_intron$Knonw[gapOl[,1]] <- bins_intron$Knonw[gapOl[,1]] - width(gap.int)[gapOl[,2]]

# get repeat info
rep = rep_info(spec1=spec1,genome=genome)
# sort into intergenic bins and intronic bins
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep), repType = rep("repeats",length(rep)))
intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep), repType = rep("repeats",length(rep)))

intergenic_chromatin <- binSort(rep = `H1-hESC`, bins = bins_gene_gap, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))
intron_chromatin <- binSort(rep = `H1-hESC`, bins = bins_intron, TE.names = names(`H1-hESC`), repType = rep("chromatin",length(`H1-hESC`)))




#### sort it out for intergenic regions
#### for each region we look for the postion of the closest 3' and 5' ends 
#### or should we use our old covCalc to construct a table for each repeat familiy 


# or maybe do it with bin sort.

# this script is still in progress



lenChoice <- 500000
repChoice <- "L1PA"
repBins <- intergenic_reps$counts
type = "intergenic"
repType = "repeats"
repList = rep
polyL1 <- polyL1coor


chromList <- `H1-hESC`



covImage <- function(lenChoice, repBins , repList ,chromList, 
                                    refgene , type ,polyL1 = NULL , minRepCov = NULL, maxRepCov = NULL, 
                                    minRepSize = NULL, maxRepSize = NULL, minBinSize = NULL, 
                                    maxBinSize = NULL){
  
  chromR <- chromList$R
  colnames(chromR)[2:4] <- c("genoName", "genoStart", "genoEnd")
  repList = c(repList, list(represedChromatin = chromR))
  if(length(polyL1) != 0){
    repList = c(repList, list(polyL1 = polyL1))
  }
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
  
  repBins <- repBins[, c("chr", "start", "end", "Known", repChoice)]
  
  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  colnames(repBins)[5] = "repCov"
  if(!is.null(maxRepCov)){
    repBins <- repBins[repBins$repCov < maxRepCov,]
  }
  if(!is.null(minRepCov)){
    repBins <- repBins[repBins$repCov > minRepCov,]
  }
  
  if(!is.null(maxRepSize)){
    repGR <- repGR[width(repGR) < maxRepSize,]
  }
  if(!is.null(minRepSize)){
    repGR <- repGR[width(repGR) > minRepSize,]
  }
  seqLen <- lenChoice + 1
  names(seqLen) <- "seq"
  
  # it might be easier if we break down the four groups at the start
  
  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))
  
  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  
  if(type == "intergenic"){
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some intergenic regions have no upstream gene")
    }
    us.strand <- data.frame(intergenicID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intergenic regions")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, maxgap=1))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some intergenic regions have no downstream gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intergenic regions")
    }
    
  }else if(type == "intron"){
    
    us.s.ol <- as.matrix(findOverlaps(us.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(us.s.ol[,1])) == length(us.bins.gr))){
      stop("some upstream intron regions don't belong to a gene")
    }
    us.strand <- data.frame(intronID = us.s.ol[,1],strand = refgene[us.s.ol[,2], 4])
    us.plus <- unique(us.strand$intronID[us.strand$strand == "+"])
    us.minus <- unique(us.strand$intronID[us.strand$strand == "-"])
    
    if(length(c(us.plus[us.plus %in% us.minus],us.minus[us.minus %in% us.plus])) != 0){
      stop("strand confusion in upstream intron region")
    }
    
    ds.s.ol <- as.matrix(findOverlaps(ds.bins.gr, refgene.gr, minoverlap=2))
    if(!(length(unique(ds.s.ol[,1])) == length(ds.bins.gr))){
      stop("some downstream intron regions don't belong to a gene")
    }
    ds.strand <- data.frame(intergenicID = ds.s.ol[,1],strand = refgene[ds.s.ol[,2], 4])
    ds.plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
    ds.minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
    
    if(length(c(ds.plus[ds.plus %in% ds.minus],ds.minus[ds.minus %in% ds.plus])) != 0){
      stop("strand confusion in downstream intron region")
    }
  }else{
    stop("type needs to be intron or intergenic")
  }

  
  us.bins.gr_5 <- us.bins.gr[us.minus]
  binStartsUs5 <- as.data.frame(us.bins.gr_5)
  binStartsUs5 <- data.frame( binStartsUs5, matrix(NA, nrow = length(us.bins.gr_5), ncol = length(repList),dimnames = list(1:length(us.bins.gr_5), names(repList))))
  covListUs5 <- NULL
  
  
  us.bins.gr_3 <- us.bins.gr[us.plus]
  binStartsUs3 <- as.data.frame(us.bins.gr_3)
  binStartsUs3 <- data.frame( binStartsUs3, matrix(NA, nrow = length(us.bins.gr_3), ncol = length(repList),dimnames = list(1:length(us.bins.gr_3), names(repList))))
  covListUs3 <- NULL
  
  
  ds.bins.gr_5 <- ds.bins.gr[ds.plus]
  binStartsDs5 <- as.data.frame(ds.bins.gr_5)
  binStartsDs5 <- data.frame( binStartsDs5, matrix(NA, nrow = length(ds.bins.gr_5), ncol = length(repList),dimnames = list(1:length( ds.bins.gr_5), names(repList))))
  covListDs5 <- NULL
  
  ds.bins.gr_3 <- ds.bins.gr[ds.minus]
  binStartsDs3 <- as.data.frame(ds.bins.gr_3)
  binStartsDs3 <- data.frame( binStartsDs3, matrix(NA, nrow = length(ds.bins.gr_3), ncol = length(repList),dimnames = list(1:length( ds.bins.gr_3), names(repList))))
  covListDs3 <- NULL
  
  
  # set the loop up here
  for(r in 1:length(repList)){
    repGR <- GRanges(seqnames=Rle(repList[[r]]$genoName),
                     ranges = IRanges(start = repList[[r]]$genoStart, end = repList[[r]]$genoEnd -1))
    us.rep.int_5 <- intersect(repGR, us.bins.gr_5)
    us.Ol_5 <- as.matrix(GenomicRanges::findOverlaps(us.bins.gr_5, us.rep.int_5, select = "first"))
    us.Ol_5 <- data.frame(1:nrow(us.Ol_5), us.Ol_5)
    us.Ol_5 <- us.Ol_5[complete.cases(us.Ol_5),]
    binStartsUs5[us.Ol_5[,1],names(repList)[r]] <- start(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1
    
    us.Ol_5 <- as.matrix(GenomicRanges::findOverlaps(us.bins.gr_5, us.rep.int_5, select = "all"))
    us.cov_5 <- data.frame(start = start(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1, 
                           end =end(us.rep.int_5[us.Ol_5[,2]]) - start(us.bins.gr_5[us.Ol_5[,1]]) + 1,
                           binID = us.Ol_5[,1], width.rank = (rank(width(us.bins.gr_5),ties.method = "random"))[us.Ol_5[,1]], width= width(us.bins.gr_5)[us.Ol_5[,1]])
    
    covListUs5 <- c(covListUs5, list(us.cov_5))
    
    
    us.rep.int_3 <- intersect(repGR, us.bins.gr_3)
    us.Ol_3 <- as.matrix(findOverlaps(us.bins.gr_3, us.rep.int_3, select = "first"))
    us.Ol_3 <- data.frame(1:nrow(us.Ol_3), us.Ol_3)
    us.Ol_3 <- us.Ol_3[complete.cases(us.Ol_3),]
    binStartsUs3[us.Ol_3[,1],names(repList)[r]] <- start(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1
    
    
    us.Ol_3 <- as.matrix(findOverlaps(us.bins.gr_3, us.rep.int_3))
    us.cov_3 <- data.frame(start = start(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1, 
                           end =end(us.rep.int_3[us.Ol_3[,2]]) - start(us.bins.gr_3[us.Ol_3[,1]]) + 1,
                           binID = us.Ol_3[,1], width.rank = (rank(width(us.bins.gr_3),ties.method = "random"))[us.Ol_3[,1]], width= width(us.bins.gr_3)[us.Ol_3[,1]])
    
    covListUs3 <- c(covListUs3, list(us.cov_3))
    
    
    
    ds.rep.int_5 <- intersect(repGR, ds.bins.gr_5)
    ds.Ol_5 <- as.matrix(findOverlaps(ds.bins.gr_5, ds.rep.int_5, select = "first"))
    ds.Ol_5 <- data.frame(1:nrow(ds.Ol_5), ds.Ol_5)
    ds.Ol_5 <- ds.Ol_5[complete.cases(ds.Ol_5),]
    binStartsUs3[ds.Ol_5[,1],names(repList)[r]] <- (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (start(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1
    
    ds.Ol_5 <- as.matrix(findOverlaps(ds.bins.gr_5, ds.rep.int_5))
    ds.cov_5 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (start(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1 ,
                           end =  (lenChoice - (end(ds.bins.gr_5[ds.Ol_5[,1]]) - start(ds.bins.gr_5[ds.Ol_5[,1]]))) + (end(ds.rep.int_5[ds.Ol_5[,2]]) - start(ds.bins.gr_5[ds.Ol_5[,1]])) + 1,
                           binID = ds.Ol_5[,1], width.rank = (rank(width(ds.bins.gr_5),ties.method = "random"))[ds.Ol_5[,1]], width= width(ds.bins.gr_5)[ds.Ol_5[,1]])
    
    covListDs5 <- c(covListDs5, list(ds.cov_5))
    
    ds.rep.int_3 <- intersect(repGR, ds.bins.gr_3)
    ds.Ol_3 <- as.matrix(findOverlaps(ds.bins.gr_3, ds.rep.int_3, select = "first"))
    ds.Ol_3 <- data.frame(1:nrow(ds.Ol_3), ds.Ol_3)
    ds.Ol_3 <- ds.Ol_3[complete.cases(ds.Ol_3),]
    binStartsUs3[ds.Ol_3[,1],names(repList)[r]] <- (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (start(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1
    
    ds.Ol_3 <- as.matrix(findOverlaps(ds.bins.gr_3, ds.rep.int_3))
    ds.cov_3 <- data.frame(start =   (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (start(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1 ,
                           end =  (lenChoice - (end(ds.bins.gr_3[ds.Ol_3[,1]]) - start(ds.bins.gr_3[ds.Ol_3[,1]]))) + (end(ds.rep.int_3[ds.Ol_3[,2]]) - start(ds.bins.gr_3[ds.Ol_3[,1]])) + 1,
                           binID = ds.Ol_3[,1], width.rank = (rank(width(ds.bins.gr_3),ties.method = "random"))[ds.Ol_3[,1]], width= width(ds.bins.gr_3)[ds.Ol_3[,1]])
    
    covListDs3 <- c(covListDs3, list(ds.cov_3))
    
    
  }
  names(covListUs5) <- names(repList)
  names(covListUs3) <- names(repList)
  names(covListDs5) <- names(repList)
  names(covListDs3) <- names(repList)
  
  
  # we need to trun some of these around to match their relative gene position
  
  
  output = list(covListUs5 = covListUs5, 
                covListDs5 = covListDs5,
                covListUs3 = covListUs3,
                covListDs3 = covListDs3,
                startsUs5 = binStartsUs5,
                startsDs5 = binStartsDs5,
                startsUs3 = binStartsUs3,
                startsDs3 = binStartsDs3
  )
  return(output)
  
}  





intronAnalysis <- covImage(lenChoice = 100000, repList = rep, chromList = `H1-hESC`, repBins = intron_reps$counts, refgene = refgene, type = "intron", polyL1 = polyL1coor)




covSamp <- intronAnalysis$covListUs5$represedChromatin
covSampStart <- intronAnalysis$startsUs5
pull = sample(1:nrow(covSampStart), 5000, replace = FALSE)
covSamp <- covSamp[covSamp$binID %in% pull, ]

width.df <- data.frame(unique(covSamp$width) ,rank(unique(covSamp$width), ties.method = "random"))

plot(c(0,100000), c(0,nrow(covSamp)), type = "n")
for(i in 1:nrow(covSamp)){
  lines(c(covSamp$start[i], covSamp$end[i]), rep(width.df[covSamp$width[i] ==width.df[,1] ,2],2), lwd = .1)
}




# how many regions of each size 
hist( log10( (chromR$genoEnd - chromR$genoStart)[(chromR$genoEnd - chromR$genoStart) %% 200 == 0] /200  ), breaks= 1000)

# heatmap the coverage 







