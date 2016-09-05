# build a table of RTNs 



# how many introns do we have 

#introns 
# introns not overlaping genes
# introns not overlaoing exons


# intergenic regions
# intergenic regions not adjacent douple stranded boundaries

filterIntron <- function(refgene){
  
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
  Echr <- NULL
  Ichr <- NULL
  Estart <- strsplit(x=as.character(refgene[,10]),split=",")
  Eend <- strsplit(x=as.character(refgene[,11]),split=",")
  Istart <- Eend
  Iend <- Estart
  for(i in 1:nrow(refgene)){
    if(refgene[i,9] > 1){
      Istart[[i]] = as.numeric(Istart[[i]][1:(length(Istart[[i]])-1)])
      Iend[[i]] = as.numeric(Iend[[i]][2:(length(Iend[[i]]))])
    }else{
      Istart[[i]] = NA
      Iend[[i]] = NA
    }
    Estart[[i]] = as.numeric(Estart[[i]])
    Eend[[i]] = as.numeric(Eend[[i]])
    Echr <- c(Echr, list(rep(as.character(refgene[i,3]), refgene[i,9])))
    Ichr <- c(Ichr, list(rep(as.character(refgene[i,3]), length(Istart[[i]]))))
  }
  
  
  
  
  Exons <- data.frame(chr = do.call(c, Echr), start = do.call(c,Estart), end = do.call(c,Eend))
  Introns <- data.frame(chr = do.call(c, Ichr), start = do.call(c,Istart), end = do.call(c,Iend))
  Introns <- Introns[!(is.na(Introns$start)),]                                        
  
  
  # intron filtering, probably could be a function
  
  intron.gr <- GRanges(seqnames = Rle(Introns$chr), 
                       ranges = IRanges(start=Introns$start, end=Introns$end)
  )
  print("number of introns")
  print(length(intron.gr))
  
  I.OL <- as.matrix(findOverlaps(intron.gr))
  intronPull <- (1:length(intron.gr))[-unique(I.OL[duplicated(I.OL[,1]),1])]
  # keep unique ones
  intronKeep.gr <- intron.gr[intronPull]
  
  intronReamin.gr <- intron.gr[unique(I.OL[duplicated(I.OL[,1]),1])]
  I.OL <- as.matrix(findOverlaps(intronReamin.gr, type = "equal"))
  eqOl.gr <- intronReamin.gr[unique(I.OL[duplicated(I.OL[,1]),1])]
  red <- reduce(eqOl.gr)
  ol2 <- as.matrix(findOverlaps(red, intronReamin.gr[-unique(I.OL[duplicated(I.OL[,1]),1])]))
  intronKeep.gr = c(intronKeep.gr, red[-unique(ol2[,1])])
  
  if(length(intronKeep.gr) == length(reduce(intronKeep.gr))){
    print("Single layer of non alt splice intronic regions")
    print(length(intronKeep.gr))
  }
  
  
  exon.gr <- GRanges(seqnames=Rle(Exons$chr),
                     ranges = IRanges(start=Exons$start, end = Exons$end))
  
  IolE <- as.matrix(findOverlaps(intronKeep.gr,exon.gr,minoverlap=2))
  
  intronKeep.gr <- intronKeep.gr[-unique(IolE[,1])]
  print("introns overlapping exons")
  print(length(intronKeep.gr))
  # strand consistancy 
  
  IolS <- as.matrix(findOverlaps(intronKeep.gr,refgene.gr,minoverlap=2))
  strandCon <- data.frame(intronID = IolS[,1],strand = refgene[IolS[,2], 4])
  plus <- unique(strandCon$intronID[strandCon$strand == "+"])
  minus <- unique(strandCon$intronID[strandCon$strand == "-"])
  strandConflict <- unique(c(plus[plus %in% minus] , minus[minus %in% plus]))
  if(length(strandConflict) > 0 ){
    intronKeep.gr <- intronKeep.gr[-strandConflict]
  }
  if(length(grep("_", seqnames(intronKeep.gr))) > 0){
    intronKeep.gr <- intronKeep.gr[-(grep("_", seqnames(intronKeep.gr)))]
  }
  return(intronKeep.gr[width(intronKeep.gr)>200])
  
}






filterIntergenic <- function(refgene){
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  refgene_gap.gr <- gaps(refgene.gr)
  refgene_gap.gr <- refgene_gap.gr[-(grep("_", seqnames(refgene_gap.gr)))]
  
  print("total intergenic regions")
  print(length(refgene_gap.gr))
  
  refgene_gap_startCH.gr <- refgene_gap.gr 
  end(refgene_gap_startCH.gr) <- end(refgene_gap.gr) - 1
  Sol <- as.matrix(findOverlaps(refgene_gap_startCH.gr, refgene.gr, maxgap=1))
  
  us.strand <- data.frame(intergenicID = Sol[,1],strand = refgene[Sol[,2], 4])
  plus <- unique(us.strand$intergenicID[us.strand$strand == "+"])
  minus <- unique(us.strand$intergenicID[us.strand$strand == "-"])
  strandConflict <- unique(c(plus[plus %in% minus] , minus[minus %in% plus]))
  
  refgene_gap_endCH.gr <- refgene_gap.gr 
  start(refgene_gap_endCH.gr) <- start(refgene_gap.gr) + 1
  Eol <- as.matrix(findOverlaps(refgene_gap_endCH.gr, refgene.gr, maxgap=1))
  
  ds.strand <- data.frame(intergenicID = Eol[,1],strand = refgene[Eol[,2], 4])
  plus <- unique(ds.strand$intergenicID[ds.strand$strand == "+"])
  minus <- unique(ds.strand$intergenicID[ds.strand$strand == "-"])
  strandConflict<- c(strandConflict, unique(c(plus[plus %in% minus] , minus[minus %in% plus])))
  
  gapKeep <- unique(c(Eol[Eol[,1] %in% Sol[,1],1],   Sol[Sol[,1] %in% Eol[,1],1]))
  print("genes on both ends")
  print(length(gapKeep))
  
  gapKeep <- gapKeep[!(gapKeep %in% strandConflict)]
  
  refgene_gap.gr <- refgene_gap.gr[gapKeep]
  #refgene_gap.gr <- refgene_gap.gr[-(grep("_", seqnames(refgene_gap.gr)))]
  print("final intergenic regions")
  print(length(refgene_gap.gr[width(refgene_gap.gr)>200]))
  return(refgene_gap.gr[width(refgene_gap.gr)>200])
}




genome = "hg19"


library(GenomicRanges)

web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)


filterIntron(refgene = refgene)

filterIntergenic(refgene)




