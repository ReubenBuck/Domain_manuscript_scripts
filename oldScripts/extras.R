# lets set up a seperate GC function 


# this one we can use to get the GC rates or at least the sequences
# then we pull the ones we want and also sample for the regression analysis 


# I'm inclined to think our GC estimates are wrong because we are 
# taking a rate of a whole interval not a section of the genome we are interested in 



# 



#lenChoice = 100000
repBins = intron_reps$counts
type = "intron"




library(BSgenome)
bsGenome <- available.genomes()[grep(genome,available.genomes())]
bsGenome <- bsGenome[-(grep("masked", bsGenome))]
library(bsGenome, character.only=TRUE)
bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])





refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))

repBins <- repBins[, c("chr", "start", "end", "Known")]






us.ends <- as.integer(repBins$start + ((repBins$end - repBins$start)/2))
us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                      ranges = IRanges(start=us.bins$start, end = us.bins$end))




ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
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




# how do we keep the order


# we could set up the matrix with the correct row names corresponding to each intron id 



All.seq<- DNAStringSet()
prime5grDS <- ds.bins[ds.plus,]
prime5grUS <- us.bins[us.minus,]

chromos <- unique(c( as.character(prime5grDS$chr),  as.character(prime5grUS$chr)))


for(i in 1:length(chromos)){
  Seq.setDS=DNAStringSet(bsSpec[[chromos[i]]], 
                         start=prime5grDS$start[prime5grDS$chr == chromos[i]], 
                         end=prime5grDS$end[prime5grDS$chr == chromos[i]])
  
  
  Seq.setUS=DNAStringSet(bsSpec[[chromos[i]]], 
                         start=prime5grUS$start[prime5grUS$chr == chromos[i]], 
                         end=as.integer(prime5grUS$end[prime5grUS$chr == chromos[i]]))
  
  All.seq <- c(All.seq, reverse(Seq.setDS), Seq.setUS)
  print(paste("5 prime", chromos[i]))
}

gPos <- vmatchPattern("C", All.seq)
#G.cov <- coverage(gPos)
G.start <- start(gPos)



Gmat <- matrix(data=NA, nrow=length(G.start), ncol=max(width(All.seq)))






gcPos <- c((vmatchPattern("C", All.seq)),(vmatchPattern("G", All.seq)))
prime5gcCov <- as.numeric(coverage(gcPos))






All.seq<- DNAStringSet()
prime3grDS <- ds.bins[ds.minus,]
prime3grUS <- us.bins[us.plus,]

chromos <- as.character(unique(c(prime3grDS$chr, prime3grUS$chr))
for(i in 1:length(chromos)){
  Seq.setDS=DNAStringSet(Hsapiens[[chromos[i]]], 
                         start=prime3grDS$start[prime3grDS$chr == chromos[i]], 
                         end=prime3grDS$end[prime3grDS$chr == chromos[i]])
  
  
  Seq.setUS=DNAStringSet(Hsapiens[[chromos[i]]], 
                         start=prime3grUS$start[prime3grUS$chr == chromos[i]], 
                         end=as.integer(prime3grUS$end[prime3grUS$chr == chromos[i]]))
  
  All.seq <- c(All.seq, reverse(Seq.setDS), Seq.setUS)
  print(paste("3 prime", chromos[i]))
}

gcPos <- c(unlist(vmatchPattern("C", All.seq)),unlist(vmatchPattern("G", All.seq)))
prime3gcCov <- as.numeric(coverage(gcPos))






# maybe we should try looking at this again with the heatmap idea 
# plot out where things tend to end up 
# it can maybe point at a length bias. 


lenChoice = 40000
repChoice = "L1PA"
repBins = intron_reps$counts[sample(x=1:nrow(intron_reps$counts), size=1000,replace=F),]
repList = rep
refgene = refgene
type = "intron"





covCalcPlot5prime3prime <- function(lenChoice, repChoice, repBins , repList , 
                                    refgene , type , minRepCov = NULL, maxRepCov = NULL, 
                                    minRepSize = NULL, maxRepSize = NULL, minBinSize = NULL, 
                                    maxBinSize = NULL){
  
  
  
  refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
  
# repBins <- repBins[repBins[,repChoice] > 0, c("chr", "start", "end", "Known", repChoice)]
repBins <- repBins[, c("chr", "start", "end", "Known", repChoice)]  


  if(!is.null(maxBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 < maxBinSize,]
  }
  if(!is.null(minBinSize)){
    repBins <- repBins[repBins$end - repBins$start + 1 > minBinSize,]
  }
  repBins2 <- repBins
  
  colnames(repBins)[5] = "repCov"
  if(!is.null(maxRepCov)){
    repBins <- repBins[repBins$repCov < maxRepCov,]
  }
  if(!is.null(minRepCov)){
    repBins <- repBins[repBins$repCov > minRepCov,]
  }
  
  repGR <- GRanges(seqnames=Rle(repList[[repChoice]]$genoName),
                   ranges = IRanges(start=repList[[repChoice]]$genoStart, end=repList[[repChoice]]$genoEnd)
  )
  
  if(!is.null(maxRepSize)){
    repGR <- repGR[width(repGR) < maxRepSize,]
  }
  if(!is.null(minRepSize)){
    repGR <- repGR[width(repGR) > minRepSize,]
  }
  seqLen <- lenChoice
  names(seqLen) <- "seq"
  
  us.ends <- repBins$start + ((repBins$end - repBins$start)/2)
  us.ends[us.ends - repBins$start+1 > lenChoice] <- repBins$start[us.ends - repBins$start+1 > lenChoice] + lenChoice
  us.bins <- data.frame(chr = repBins$chr, start = repBins$start, end = us.ends)
  us.bins.gr <- GRanges(seqnames=Rle(us.bins$chr), 
                        ranges = IRanges(start=us.bins$start, end = us.bins$end))
  us.rep.int <- intersect(repGR, us.bins.gr)
  us.Ol <- as.matrix(findOverlaps(us.bins.gr, us.rep.int))
  us.cov <- data.frame(start = start(us.rep.int[us.Ol[,2]]) - start(us.bins.gr[us.Ol[,1]]) + 1, 
                       end =end(us.rep.int[us.Ol[,2]]) - start(us.bins.gr[us.Ol[,1]]) + 1)
  
  
  
  ds.starts <- as.integer(repBins$end - ((repBins$end - repBins$start)/2))
  ds.starts[repBins$end - ds.starts + 1> lenChoice] <- repBins$end[repBins$end - ds.starts + 1> lenChoice] - lenChoice
  ds.bins <- data.frame(chr = repBins$chr, start = ds.starts, end = repBins$end)
  ds.bins.gr <- GRanges(seqnames=Rle(ds.bins$chr), 
                        ranges = IRanges(start=ds.bins$start, end = ds.bins$end))
  ds.rep.int <- intersect(repGR, ds.bins.gr)  
  ds.Ol <- as.matrix(findOverlaps(ds.bins.gr, ds.rep.int))
  ds.cov <- data.frame(start =   (lenChoice - (end(ds.bins.gr[ds.Ol[,1]]) - start(ds.bins.gr[ds.Ol[,1]]))) + (start(ds.rep.int[ds.Ol[,2]]) - start(ds.bins.gr[ds.Ol[,1]])) + 1 ,
                       end =  (lenChoice - (end(ds.bins.gr[ds.Ol[,1]]) - start(ds.bins.gr[ds.Ol[,1]]))) + (end(ds.rep.int[ds.Ol[,2]]) - start(ds.bins.gr[ds.Ol[,1]])) + 1)
  
  
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
  ### everything that is minus we have to reverse the coordinates 
  # upstream region of intergenic in plus strand is the 3prime region of a gene
  #us.Ol[,1] tells us which bin the repeats belong to
  prime3_us <- us.cov[us.Ol[us.Ol[,1] %in% us.plus,2],]
  prime5_us <- data.frame(start = lenChoice + 1 - us.cov[us.Ol[us.Ol[,1] %in% us.minus,2],]$end + 1, 
                          end = lenChoice + 1 - us.cov[us.Ol[us.Ol[,1] %in% us.minus,2],]$start + 1)
  
  prime3_ds <- data.frame(start = lenChoice + 1- ds.cov[ds.Ol[ds.Ol[,1] %in% ds.minus,2],]$end + 1, 
                          end = lenChoice + 1 - ds.cov[ds.Ol[ds.Ol[,1] %in% ds.minus,2],]$start + 1)
  prime5_ds <- ds.cov[ds.Ol[ds.Ol[,1] %in% ds.plus,2],]
  
  prime5 <- rbind(prime5_ds, prime5_us)
  prime3 <- rbind(prime3_ds, prime3_us)
  
  prime5.cov.r <- GRanges(seqnames=Rle("seq"), ranges = IRanges(start=prime5$start, end = prime5$end), seqlengths = seqLen + 1)
  prime3.cov.r <- GRanges(seqnames=Rle("seq"), ranges = IRanges(start=prime3$start, end = prime3$end), seqlengths = seqLen + 1)
  
  #### coverage is sorted. 
  ### now i need to get base frequencies for each side. 
  ### this means collecting the correct information 
  ### us.plus and us.minus can help sort that
  
  
  prime5Lengths <- c(width(us.bins.gr[us.minus]), width(ds.bins.gr[ds.plus]))
  prime3Lengths <- c(width(us.bins.gr[us.plus]), width(ds.bins.gr[ds.minus]))
  
  len.s <- sort(prime5Lengths)
  lenSum <- summary(as.factor(len.s), maxsum=length(unique(len.s)))
  CS <- cumsum(as.integer(lenSum)[length(lenSum):1])
  CS <- CS[length(CS):1]
  step <- stepfun(as.numeric(names(lenSum)[2:(length(lenSum))]), c(CS))
  prime5stepRes <- step(1:(lenChoice+1))[(lenChoice+1):1]
  
  len.s <- sort(prime3Lengths)
  lenSum <- summary(as.factor(len.s), maxsum=length(unique(len.s)))
  CS <- cumsum(as.integer(lenSum)[length(lenSum):1])
  CS <- CS[length(CS):1]
  step <- stepfun(as.numeric(names(lenSum)[2:(length(lenSum))]), c(CS))
  prime3stepRes <- step(1:(lenChoice+1))
  
  
  p = sum(repBins$repCov)/sum(repBins$end - repBins$start + 1)
  output = list(rawRepCov5 = as.numeric(coverage(prime5.cov.r)$seq),
                rawRepCov3 = as.numeric(coverage(prime3.cov.r)$seq) ,
                zRepCov5 = (as.numeric(coverage(prime5.cov.r)$seq) - prime5stepRes*p)/sqrt(prime5stepRes*p*(1-p)),
                zRepCov3 = (as.numeric(coverage(prime3.cov.r)$seq) - prime3stepRes*p)/sqrt(prime3stepRes*p*(1-p)), 
                p = p, 
                baseFreq5prime = prime5stepRes, 
                baseFreq3prime = prime3stepRes
  )
  return(output)
  
}  







# lets pick say 1000 random bins
# we'll do this on the upstream 

#binSamp <- sample(1:nrow(us.bins), size=100,replace=FALSE)
#cov.samp <- us.cov[us.Ol[us.Ol[,1] %in% binSamp,2]]


head(us.cov)
head(us.Ol)
head(us.bins)

mat <- matrix(1, nrow = nrow(us.bins), ncol = 40000 + 2)

for(i in 1:nrow(mat)){
  mat[i,1:(us.bins$end[i] - us.bins$start[i] + 1)] <- 0
  
  
  coords <- us.cov[us.Ol[,1] == i,]
  
  
  if(nrow(coords) > 0){
    for(c in 1:nrow(coords)){
      mat[i,coords$start[c]: coords$end[c]] <- 2
    }
  }
}



image(t(mat[order((us.bins$end - us.bins$start + 1)),]), col = grey(level=c(0,.5,1)))



118277619 - 118274985







########## Motif contnent

#### This all links up with the fig2 functions after plotting out all the raw and modified counts



# we known TE content changes as a result of GC,
# but also GC changes as a result of hterochromatin 
# if heterochromatin adjustments are able to stabalise GC it works out well

# so now we can normalize for this 
# two step normalization, this should reval what is happening



# if we balance the GC content based on the chromatin state 
# then we can balance the TE content on the new GC. 



bigPca <- prcomp(joinSampGenomeBig[,5:ncol(joinSampGenomeBig)], scale. = T)
reuben.biplot(bigPca$x,bigPca$rotation,text.col = 2,cex = .5)
smallPca <- prcomp(joinSampGenomeSmall[,9:ncol(joinSampGenomeSmall)], scale. = F)
reuben.biplot(smallPca$x,smallPca$rotation,text.col = 2,cex = .5)
L1_pca <- prcomp(joinSampGenomeBig[,c("new_L1", "GC", "R")], scale.=T,center = F)
reuben.biplot(L1_pca$x[,c(1,3)], L1_pca$rotation[,c(1,3)], text.col = 2,cex = .5)

bigMotif <- LocalMotifLevel(repBins = joinSampGenomeBig,genome = genome, motif = "TTTTAA")
smallMotif <- LocalMotifLevel(repBins = joinSampGenomeSmall,genome = genome, motif = "TTTTAA")
medMotif <- LocalMotifLevel(repBins = joinSampGenomeMed,genome = genome, motif = "TTTTAA")
largeMotif <- LocalMotifLevel(repBins = joinSampGenomeLarge,genome = genome, motif = "TTTTAA")


y = bigMotif$motif/20000
x = joinSampGenomeBig$GC*100
plot(x,y,  pch = 16, cex = .5)
lmPolyBig <- lm(y ~ x+ I(x^2))
lines(seq(0,100,1), predict(object = lmPolyBig, data.frame(x = seq(0,100,1))), col = "red",lwd = 2)

y = largeMotif$motif/30000
x = joinSampGenomeLarge$GC*100
plot(x,y,  pch = 16, cex = .5)
lmPolyLarge <- lm(y ~ x+ I(x^2))
lines(seq(0,100,1), predict(object = lmPolyLarge, data.frame(x = seq(0,100,1))), col = "red",lwd = 2)


y = medMotif$motif/10000
x = joinSampGenomeMed$GC*100
plot(x,y,  pch = 16, cex = .5)
lmPolyMed <- lm(y ~ x+ I(x^2))
lines(seq(0,100,1), predict(object = lmPolyMed, data.frame(x = seq(0,100,1))), col = "red",lwd = 2)


y = smallMotif$motif/1000
x = joinSampGenomeSmall$GC*100
plot(x,y,  pch = 16, cex = .5)
lmPolySmall <- lm(y ~ x + I(x^2))
lines(seq(0,100,1), predict(object = lmPolySmall, data.frame(x = seq(0,100,1))), col = "red",lwd = 2)


plot(1,1,ylim = c(0,.006), xlim = c(30,60), type = "n", xlab = "GC%", ylab = "insertion motif fraction")
y = c(predict(object = lmPolySmall, data.frame(x = seq(0,100,1))) - summary.lm(lmPolySmall)$sigma,
      predict(object = lmPolySmall, data.frame(x = seq(0,100,1)))[101:1] + summary.lm(lmPolySmall)$sigma)
x = c(seq(0,100,1), seq(0,100,1)[101:1])
polygon(x=x,y=y, col = "red", density = 30)
lines(seq(0,100,1), predict(object = lmPolySmall, data.frame(x = seq(0,100,1))), col = "red",lwd = 3)
y = c(predict(object = lmPolyMed, data.frame(x = seq(0,100,1))) - summary.lm(lmPolyMed)$sigma,
      predict(object = lmPolyMed, data.frame(x = seq(0,100,1)))[101:1] + summary.lm(lmPolyMed)$sigma)
x = c(seq(0,100,1), seq(0,100,1)[101:1])
polygon(x=x,y=y, col = "darkgreen", density = 30)
lines(seq(0,100,1), predict(object = lmPolyMed, data.frame(x = seq(0,100,1))), col = "darkgreen",lwd = 3)
y = c(predict(object = lmPolyBig, data.frame(x = seq(0,100,1))) - summary.lm(lmPolyBig)$sigma,
      predict(object = lmPolyBig, data.frame(x = seq(0,100,1)))[101:1] + summary.lm(lmPolyBig)$sigma)
x = c(seq(0,100,1), seq(0,100,1)[101:1])
polygon(x=x,y=y, col = "blue", density = 30)
lines(seq(0,100,1), predict(object = lmPolyLarge, data.frame(x = seq(0,100,1))), col = "blue",lwd = 3)
y = c(predict(object = lmPolyLarge, data.frame(x = seq(0,100,1))) - summary.lm(lmPolyLarge)$sigma,
      predict(object = lmPolyLarge, data.frame(x = seq(0,100,1)))[101:1] + summary.lm(lmPolyLarge)$sigma)
x = c(seq(0,100,1), seq(0,100,1)[101:1])
polygon(x=x,y=y, col = "orange", density = 30)
lines(seq(0,100,1), predict(object = lmPolyLarge, data.frame(x = seq(0,100,1))), col = "orange",lwd = 3)

legend("topright", legend = c(1000,10000,20000,30000), fill = c("red", "darkgreen", "blue", "orange"),title = "bin size")





######## more extras






heteroBig <- localGCgenome(repList = `H1-hESC`,binSize = 20000,sampSize = 5000, genome = genome, repType = "chromatin")

layout(1)
par(mar=c(5,5,5,5))
plot(heteroBig$GC, heteroBig$T, pch = 16, cex = .3)
heterolo <- loess(heteroBig$T~heteroBig$GC)
lines(seq(.1,1,.01), predict(heterolo, seq(.1,1,.01)))


Rlo  <- loess(heteroBig$R~heteroBig$GC)
Elo  <- loess(heteroBig$E~heteroBig$GC)
PFlo <- loess(heteroBig$PF~heteroBig$GC)
CTCFlo <- loess(heteroBig$CTCF~heteroBig$GC)
Ttlo <- loess(heteroBig$T~heteroBig$GC)
TSSlo <- loess(heteroBig$TSS~heteroBig$GC)
WElo <- loess(heteroBig$WE~heteroBig$GC)



# we need to get bins that have chromatin status included 


colnames(intergenic_chromatin$counts)[5:ncol(intergenic_chromatin$counts)]

lenChoice = 60000

R <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="R", repBins=intergenic_chromatin$counts[intergenicSample,],repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")

Rsum5 <- rollsum(R$rawRepCov5,k = k, fill="extend", align = "center")
Rrate <- Rsum5/baseSum5

gc5adj <- gcRate5 - predict(Rlo, Rrate)
plot(gc5adj + mean(joinSampGenomeBig$GC), type = "l", ylim = c(0.3,0.6))

lines(gcRate5,col = 2)



E <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="E", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
PF <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="PF", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
CTCF <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="CTCF", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
Tt <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="T", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
TSS <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="TSS", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")
WE <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice="WE", repBins=intergenic_chromatin$counts,repList=`H1-hESC`, refgene=refgene, type= genome_type, repType = "chromatin")




layout(matrix(c(1,2), nrow=1))

par(mar=c(5,5,5,0))
plot(-(sqrt(100001:1)),((R$rawRepCov5/R$baseFreq5prime) - predict(Rlo, gcRate5)) + R$p, type = "l", ylim = c(0,1), ylab = "",xaxt="n", xlab = "")
lines(-(sqrt(100001:1)),((E$rawRepCov5/E$baseFreq5prime)- predict(Elo, gcRate5)) + E$p, col = 2)
lines(-(sqrt(100001:1)),((PF$rawRepCov5/PF$baseFreq5prime)- predict(PFlo, gcRate5)) + PF$p, col = 3)
lines(-(sqrt(100001:1)),((CTCF$rawRepCov5/CTCF$baseFreq5prime)- predict(CTCFlo, gcRate5)) + CTCF$p, col = 4)
lines(-(sqrt(100001:1)),((Tt$rawRepCov5/Tt$baseFreq5prime)- predict(Ttlo, gcRate5)) + Tt$p, col = 5)
lines(-(sqrt(100001:1)),((TSS$rawRepCov5/TSS$baseFreq5prime)- predict(TSSlo, gcRate5)) + TSS$p, col = 6)
lines(-(sqrt(100001:1)),((WE$rawRepCov5/WE$baseFreq5prime)- predict(WElo, gcRate5)) + WE$p, col = 7)
par(new=TRUE)
new_L1_Z <- (new_L1$rawRepCov5/new_L1$baseFreq5prime - predict(new_L1Mod, gcRate5))   / ((1 * sqrt(n5 * predict(new_L1Mod, gcRate5) *  (1 - predict(new_L1Mod, gcRate5))))/n5)
plot(-(sqrt(100001:1)), new_L1_Z, ylim = c(-10,10) , col = 8, type = "l", axes = FALSE, ylab = "", xlab = "")
abline(h=2, lty=2)
abline(h=-2, lty=2)

axis(1,at = seq(0,-350, by = -50), labels = seq(0,350, by = 50)^2)



par(mar=c(5,0,5,5))

plot((sqrt(1:100001)),((R$rawRepCov3/R$baseFreq3prime)- predict(Rlo, gcRate3)) + R$p, type = "l", ylim = c(0,1), xaxt = "n", yaxt="n", ylab = "n", xlab = "")
lines((sqrt(1:100001)),((E$rawRepCov3/E$baseFreq3prime)- predict(Elo, gcRate3)) + E$p, col = 2)
lines((sqrt(1:100001)),((PF$rawRepCov3/PF$baseFreq3prime)- predict(PFlo, gcRate3)) + PF$p, col = 3)
lines((sqrt(1:100001)),((CTCF$rawRepCov3/CTCF$baseFreq3prime)- predict(CTCFlo, gcRate3)) + CTCF$p, col = 4)
lines((sqrt(1:100001)),((Tt$rawRepCov3/Tt$baseFreq3prime)- predict(Ttlo, gcRate3)) + Tt$p, col = 5)
lines((sqrt(1:100001)),((TSS$rawRepCov3/TSS$baseFreq3prime)- predict(TSSlo, gcRate3)) + TSS$p, col = 6)
lines((sqrt(1:100001)),((WE$rawRepCov3/WE$baseFreq3prime)- predict(WElo, gcRate3)) + WE$p, col = 7)
par(new=TRUE)
new_L1_Z <- (new_L1$rawRepCov3/new_L1$baseFreq3prime - predict(new_L1Mod, gcRate3))   / ((1 * sqrt(n3 * predict(new_L1Mod, gcRate3) *  (1 - predict(new_L1Mod, gcRate3))))/n3)
plot((sqrt(1:100001)),(new_L1_Z), ylim = c(-10,10)  , col = 8, type = "l",xlab = "", axes=FALSE, ylab = "")
abline(h=2, lty=2)
abline(h=-2, lty=2)

axis(1,at = seq(0,350, by = 50), labels = seq(0,350, by = 50)^2)
axis(4,at=seq(-10,10), labels=seq(-10,10))





# how can we have TEs with less than the expected TE rate the whole time 
## so need to proberly normalize for GC content
## lets get the rolling windo wto work properly










# region comparison

minS <- NULL
maxS <- NULL
minRepCov <- 5000
maxRepCov <- 7000
lenChoice = 100000
repChoice <- "L1PA"
TEintergenic <- covCalcPlot5prime3prime(repType = "repeats",lenChoice=lenChoice,repChoice=repChoice, repBins=intergenic_reps$counts,repList=rep, refgene=refgene, type= "intergenic",minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
TEintron <- covCalcPlot5prime3prime(repType = "repeats",lenChoice=lenChoice,repChoice=repChoice, repBins=intron_reps$counts,repList=rep, refgene=refgene, type= "intron", minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)

logCoordinates <- c(-(log10(1:length(TEintergenic$zRepCov5))[length(TEintergenic$zRepCov5):1]), 0, log10(1:length(TEintergenic$zRepCov3)))
Coordinates <- -(length(TEintergenic$zRepCov5)):length(TEintergenic$zRepCov3)

plot(logCoordinates, c(TEintergenic$zRepCov5,NA ,TEintergenic$zRepCov3), type = "l", 
     ylim = c(min(c(TEintergenic$zRepCov5,TEintergenic$zRepCov3, 
                    TEintron$zRepCov5, TEintron$zRepCov3)),
              max(c(TEintergenic$zRepCov5,TEintergenic$zRepCov3, 
                    TEintron$zRepCov5, TEintron$zRepCov3))
     )
)
lines(logCoordinates, c(TEintron$zRepCov5,NA ,TEintron$zRepCov3), col = 2)

abline(h=2,lty=2)
abline(h=-2,lty=2)



# accumulation comparison intergenic
# we look at families fo rthe convayerbelt effct
# we could look for sites with both quite easily 


## ALU

minS <- 50000
maxS <- NULL
minRepCov <- NULL
maxRepCov <- NULL
lenChoice = 50000
genome_type <- "intergenic"
binnedGenome <- get(paste(genome_type, "_reps", sep = ""))$counts


repChoice <- "AluJ"
AluJ <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "AluS"
AluS <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "AluY"
AluY <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)


logCoordinates <- c(-(log10(1:length(AluJ$zRepCov5))[length(AluJ$zRepCov5):1]), 0, log10(1:length(AluJ$zRepCov3)))
Coordinates <- -(length(AluJ$zRepCov5)):length(AluJ$zRepCov3)

plot(logCoordinates, c(AluJ$zRepCov5,NA ,AluJ$zRepCov3), type = "l", 
     ylim = c(min(c(AluJ$zRepCov5,AluJ$zRepCov3,
                    AluS$zRepCov5,AluS$zRepCov3,
                    AluY$zRepCov5,AluY$zRepCov3)),
              max(c(AluJ$zRepCov5,AluJ$zRepCov3,
                    AluS$zRepCov5,AluS$zRepCov3,
                    AluY$zRepCov5,AluY$zRepCov3))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)
lines(logCoordinates, c(AluS$zRepCov5,NA ,AluS$zRepCov3), col = 2)
lines(logCoordinates, c(AluY$zRepCov5,NA ,AluY$zRepCov3), col = 3)
legend("bottomright", c("AluJ", "AluS", "AluY"), fill = c(1,2,3))

abline(h=2,lty=2)
abline(h=-2,lty=2)





#### L1 Families


minS <- NULL
maxS <- NULL
minRepCov <- 5000
maxRepCov <- 7000
lenChoice = 50000
genome_type <- "intergenic"
binnedGenome <- get(paste(genome_type, "_reps", sep = ""))$counts

repChoice <- "L1ME"
L1ME <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MD"
L1MD <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MC"
L1MC <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MB"
L1MB <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1MA"
L1MA <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1PB"
L1PB <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)
repChoice <- "L1PA"
L1PA <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=rep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=minRepCov, maxRepCov=maxRepCov)

logCoordinates <- c(-(log10(1:length(L1ME$zRepCov5))[length(L1ME$zRepCov5):1]), 0, log10(1:length(L1ME$zRepCov3)))
Coordinates <- -(length(L1ME$zRepCov5)):length(L1ME$zRepCov3)

plot(Coordinates, c(L1MA$zRepCov5,NA ,L1MA$zRepCov3), type = "n", 
     ylim = c(min(c(L1ME$zRepCov5,L1ME$zRepCov3,
                    L1MD$zRepCov5,L1MD$zRepCov3,
                    L1MC$zRepCov5,L1MC$zRepCov3,
                    L1MB$zRepCov5,L1MB$zRepCov3,
                    L1MA$zRepCov5,L1MA$zRepCov3,
                    L1PB$zRepCov5,L1PB$zRepCov3,
                    L1PA$zRepCov5,L1PA$zRepCov3
     )),
     max(c(L1ME$zRepCov5,L1ME$zRepCov3,
           L1MD$zRepCov5,L1MD$zRepCov3,
           L1MC$zRepCov5,L1MC$zRepCov3,
           L1MB$zRepCov5,L1MB$zRepCov3,
           L1MA$zRepCov5,L1MA$zRepCov3,
           L1PB$zRepCov5,L1PB$zRepCov3,
           L1PA$zRepCov5,L1PA$zRepCov3
     ))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)
legend("bottomright", c("L1ME", "L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA"), fill = c(1,2,3,4,5,6,8))

abline(h=2,lty=2)
abline(h=-2,lty=2)

lines(Coordinates, c(L1ME$zRepCov5,NA ,L1ME$zRepCov3), col = 1)
lines(Coordinates, c(L1MD$zRepCov5,NA ,L1MD$zRepCov3), col = 2)
lines(Coordinates, c(L1MC$zRepCov5,NA ,L1MC$zRepCov3), col = 3)
lines(Coordinates, c(L1MB$zRepCov5,NA ,L1MB$zRepCov3), col = 4)
lines(Coordinates, c(L1MA$zRepCov5,NA ,L1MA$zRepCov3), col = 5)
lines(Coordinates, c(L1PB$zRepCov5,NA ,L1PB$zRepCov3), col = 6)
lines(Coordinates, c(L1PA$zRepCov5,NA ,L1PA$zRepCov3), col = 8)









### lets try joiing families

intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                                  old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                                  new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                                  Alu = rowSums(intergenic_reps$counts[,c("AluS", "AluY", "AluJ")]),
                                  Ancient = rowSums(intergenic_reps$counts[,c("MIR", "L2")])
)

intronJoinedFam <- data.frame(intron_reps$counts[,c("chr", "start", "end","Known")],
                              old_L1 = rowSums(intron_reps$counts[,c("L1ME", "L1MD", "L1MC","L1MB")]),
                              new_L1 = rowSums(intron_reps$counts[,c("L1MA", "L1PB", "L1PA","L1HS")]),
                              Alu = rowSums(intron_reps$counts[,c("AluS", "AluY", "AluJ")]),
                              Ancient = rowSums(intron_reps$counts[,c("MIR", "L2")])
)


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))


minS <- NULL
maxS <- NULL
minRepCov <- NULL
maxRepCov <- NULL
lenChoice = 50000
genome_type <- "intron"
binnedGenome <- get(paste(genome_type , "JoinedFam", sep = ""))



repChoice <- "old_L1"
old_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "new_L1"
new_L1 <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Alu"
Alu <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)
repChoice <- "Ancient"
Ancient <- covCalcPlot5prime3prime(lenChoice=lenChoice,repChoice=repChoice, repBins=binnedGenome,repList=joinRep, refgene=refgene, type= genome_type,minBinSize=minS, maxBinSize=maxS,minRepCov=NULL, maxRepCov=NULL)

logCoordinates <- c(-(log10(1:(lenChoice + 1))[(lenChoice + 1):1]), 0, log10(1:(lenChoice + 1)))
Coordinates <- -((lenChoice + 1)):(lenChoice + 1)


plot(logCoordinates, c(old_L1$zRepCov5,NA ,old_L1$zRepCov3), type = "l", 
     ylim = c(min(c(old_L1$zRepCov5,old_L1$zRepCov3,
                    new_L1$zRepCov5,new_L1$zRepCov3,
                    Alu$zRepCov5,Alu$zRepCov3,
                    Ancient$zRepCov5,Ancient$zRepCov3)),
              max(c(old_L1$zRepCov5,old_L1$zRepCov3,
                    new_L1$zRepCov5,new_L1$zRepCov3,
                    Alu$zRepCov5,Alu$zRepCov3,
                    Ancient$zRepCov5,Ancient$zRepCov3))
     ),
     main = genome_type, ylab = "TE accumulation (Z score)"
)

lines(logCoordinates, c(new_L1$zRepCov5,NA ,new_L1$zRepCov3), col = 2)
lines(logCoordinates, c(Alu$zRepCov5,NA ,Alu$zRepCov3), col = 3)
lines(logCoordinates, c(Ancient$zRepCov5,NA ,Ancient$zRepCov3), col = 4)


legend("bottomright", names(joinRep), fill = c(1,2,3,4))

abline(h=2,lty=2)
abline(h=-2,lty=2)








plot((1:length(ALuSintron$zRepCov)), ALuSintron$zRepCov,col = 1, type = "l")
lines((1:length(L1MEintron$zRepCov)), L1MEintron$zRepCov,col = 2, type = "l")



abline(h = 2)


# maybe new insertions are accoicated with a specific type of histone mark 

plot((1:length(ALuSintron$zRepCov)), ALuSintron$rawRepCov, type = "l")
lines((1:length(ALuSintron$zRepCov)), ALuS$rawRepCov, col = 2)
par(new = T)
plot((1:length(ALuSintron$zRepCov)), ALuSintron$rawRepCov/ALuSintron$baseFreq, type = "l")
lines((1:length(ALuSintron$zRepCov)), ALuS$rawRepCov/ALuSintergenic$baseFreq, col= 2)


par(new=T)
plot((1:length(ALuSintron$zRepCov)), ALuSintron$zRepCov, type = "l", col = 1)
lines((1:length(ALuSintron$zRepCov)), ALuSintergenic$zRepCov, type = "l", col = 2)

# as far as i can tell, enrichment tends to be at these smaller sizes. 


plot((1:length(ALuSintron$zRepCov)), ALuSintron$baseFreq/sum(ALuSintron$baseFreq), type = "l")
lines((1:length(ALuSintergenic$zRepCov)), ALuSintergenic$baseFreq/sum(ALuSintergenic$baseFreq), type = "l", col = 2)


plot((1:length(ALuSintron$zRepCov)), ALuSintron$baseFreq, type = "l")
lines((1:length(ALuSintergenic$zRepCov)), ALuSintergenic$baseFreq, type = "l", col = 2)









# There may be an effect where shorter intorns tend to be from shorter genes that 
# tend to get expressed more. This means they are more chromaticly open 
# leading to accumulation of TEs and a confounding factor in our results
# One straturdy is to stratify on Intron size. 
# I should cut out small introns and then aother patterns may begin to immerge

# if we are going to do the bidirectional plot we need to incorperate strand information. 
# Essentially know if something is up or downstream of a gene 
# also this will change the number of up and downstream regions we are looking at





TEchioce = "AluY"
datType = "counts"

plot(log10(intron_reps$rates$Known),  log10(intron_reps[[datType]][,TEchioce]) , pch = 16, cex = .1, col = 2, main = TEchioce,
     ylab = paste("log10", datType), xlab = "log10 interval size")

points(log10(intergenic_reps$rates$Known),log10(intergenic_reps[[datType]][,TEchioce]),  pch = 16, cex = .1)
legend("topleft", legend = c("intergenic", "intron"), fill = c(1,2))
a = (1:100)
b = (1:100)

fit <- lm(b~a)
abline(lm(b~a))



plot(log10(intron_reps$rates$Known -(intron_reps[[datType]][,TEchioce]) )  ,  log10(intron_reps[[datType]][,TEchioce]) , pch = 16, cex = .1, col = 2, main = TEchioce,
     ylab = paste("log10", datType), xlab = "log10 interval size")
points(log10( (intergenic_reps$rates$Known) -intergenic_reps[[datType]][,TEchioce] ) ,log10(intergenic_reps[[datType]][,TEchioce]),  pch = 16, cex = .1)



# can these differences just be a result of more smaller spaces

# should we stratify on region size

# also we are thinking about similar outcomes for both regions but by different mechanisms. 


plot(density(log10(intron_reps$counts$Known)), col = 2, xlab = "log10 interval length")
lines(density(log10(intergenic_reps$counts$Known)), col = 1)
legend("topleft", legend=c("intergenic", "intron"), fill=c(1,2))


bin_strat <- c(seq(100, 900, by = 100), seq(1000, 9000, by = 1000), seq(10000, 100000, by = 10000))
result_tabInter <- data.frame(intervalSize = bin_strat[1:(length(bin_strat)-1)], mean = NA, sd = NA)
result_tabIntron <- data.frame(intervalSize = bin_strat[1:(length(bin_strat)-1)], mean = NA, sd = NA)


repChoice = "AluS"
InterDat = NULL
for(i in 1:nrow(result_tabInter)){
  selection = intergenic_reps$rates[,repChoice][intergenic_reps$rates$Known < bin_strat[i+1] & intergenic_reps$rates$Known > bin_strat[i]]
  selection = selection[selection > 0]
  result_tabInter$mean[i] = mean(selection)
  result_tabInter$sd[i] = sd(selection)
  InterDat <- c(InterDat, list(selection))
}

IntronDat = NULL
for(i in 1:nrow(result_tabIntron)){
  selection = intron_reps$rates[,repChoice][intron_reps$rates$Known < bin_strat[i+1] & intron_reps$rates$Known > bin_strat[i]]
  selection = selection[selection > 0]
  result_tabIntron$mean[i] = mean(selection)
  result_tabIntron$sd[i] = sd(selection)
  IntronDat <- c(IntronDat, list(selection))
  
}
names(InterDat) <- result_tabInter$intervalSize
names(IntronDat) <- result_tabIntron$intervalSize



plot(log10(result_tabInter$intervalSize), result_tabInter$mean, type = "b", ylim=c(1,10000))
lines(log10(result_tabInter$intervalSize), result_tabInter$mean + (result_tabInter$sd), lty = 2)
lines(log10(result_tabInter$intervalSize), result_tabInter$mean-(result_tabInter$sd), lty = 2)

lines(log10(result_tabIntron$intervalSize), result_tabIntron$mean, type = "b", col = 2)
lines(log10(result_tabIntron$intervalSize), result_tabIntron$mean + (result_tabIntron$sd), lty = 2, col = 2)
lines(log10(result_tabIntron$intervalSize),  result_tabIntron$mean - (result_tabIntron$sd), lty = 2, col = 2)


boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)
par(new = T)
boxplot(InterDat,outline=T, ylim = c(1,10000))


boxplot(InterDat,outline=T, ylim = c(1,10000))
par(new = T)
boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)


boxplot(IntronDat,outline=T, ylim = c(1,10000), col = 2)
boxplot(InterDat,outline=T, ylim = c(1,10000))




# other experiments to try will be to look at heterochromatin coverage
# and DNase1 coverage
# also coorect for GC content

# not too hard, there is a bioconductor packakage we can use to get that info. 
# we could potentially do this today 
# we know which ranges we are grabbing for both intron and intergenic 
# all we need is to factor it in
### lets try pull out GC content

Seq=Hsapiens[["chr1"]]


chrChoice <- "chr21"

Seq <- Hsapiens[["chr21"]]

# we can do this for each region via chromosome. 
Seq.set=DNAStringSet(Hsapiens[[chrChoice]], start=start(refgene_gap.gr[seqnames(refgene_gap.gr) == chrChoice]), end=end(refgene_gap.gr[seqnames(refgene_gap.gr) == chrChoice]))
bases=alphabetFrequency(Seq.set,baseOnly=TRUE)
CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])



# find 

chromos <- as.character(unique(intergenic_reps$counts$chr))
IntergenicGC <- rep(NA, nrow(intergenic_reps$counts))
for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=intergenic_reps$counts$start[intergenic_reps$counts$chr == chromos[i]], 
                       end=intergenic_reps$counts$end[intergenic_reps$counts$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  IntergenicGC[intergenic_reps$counts$chr == chromos[i]] <- CGcontent
  
}

IntronGC <- rep(NA, nrow(intron_reps$counts))
for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=intron_reps$counts$start[intron_reps$counts$chr == chromos[i]], 
                       end=intron_reps$counts$end[intron_reps$counts$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  IntronGC[intron_reps$counts$chr == chromos[i]] <- CGcontent
  
}


repChoice <- "AluS"
intronGCset <- IntronGC[intron_reps$rates[,repChoice] > 0]
intronTEset <- intron_reps$rates[intron_reps$rates[,repChoice] > 0,repChoice]
intergenicGCset <- IntergenicGC[intergenic_reps$rates[,repChoice] > 0]
intergenicTEset <- intergenic_reps$rates[intergenic_reps$rates[,repChoice] > 0,repChoice]

plot((intronGCset),
     (intronTEset),
     pch = 16, cex = .05, col = 2 , main = repChoice,
     ylab = "TE coverage rate (log10 per 10000bp)",
     xlab = "GC content"
)
points((intergenicGCset),
       (intergenicTEset), 
       pch = 16, cex = .05)
interMod <- lm((intergenicTEset) ~ (intergenicGCset))
intronMod <- lm((intronTEset) ~ (intronGCset))


intronMod2 <- loess((intronTEset) ~ (intronGCset))
interMod2 <- loess((intergenicTEset) ~ (intergenicGCset))
lines(seq(.1,.65,.01), predict(object=intronMod2, newdata=seq(.1,.65,.01)), col = 2, lwd = 2)
lines(seq(.1,.65,.01), predict(object=interMod2, newdata=seq(.1,.65,.01)), col = 1, lwd = 2)



abline(intronMod, col = 2, lwd = 2, lty = 2)
abline(interMod, col = 1, lwd = 2, , lty = 2)
legend("bottomright", legend=c("intergenic", "intron"), fill = c(1,2))

legend("topright", bty="n", legend=c(paste("intergenic r =", format(sqrt(summary(interMod)$r.squared), digits=3)), 
                                     paste("intron r =", format(sqrt(summary(intronMod)$r.squared), digits=3))
)
)




### if we can extract the Cs and Gs and put them in a rep style list opject
### then we can plug them straight into our annalysis 
### the Cs and Gs can be in any order



# this is essentially another way to work it out 
# it is probably quiclker then calculating each one individually and then sending it through our script
# alternativly 
# we could just add the start position of the interval 
# this function gcRepObject


# lets see if we can vectorise when we calculate the mean 
chromos <- as.character(unique(repBins$chr))
repBinsGC <- rep(NA, nrow(repBins))
for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=repBins$start[repBins$chr == chromos[i]], 
                       end=repBins$end[repBins$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  repBinsGC[repBins$chr == chromos[i]] <- CGcontent
  
}



chromos <- as.character(unique(repBins$chr))
repBinsGC <- rep(NA, nrow(repBins))


for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(Hsapiens[[chromos[i]]], 
                       start=repBins$start[repBins$chr == chromos[i]], 
                       end=repBins$end[repBins$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  repBinsGC[repBins$chr == chromos[i]] <- CGcontent
  
}


# calculate residual 


# maybe just get bins 10,000 bp 
# do the regression analysis for each repeat 



library(BSgenome)
bsGenome <- available.genomes()[grep(genome,available.genomes())]
bsGenome <- bsGenome[-(grep("masked", bsGenome))]
library(bsGenome, character.only=TRUE)
bsSpec <- get(strsplit(bsGenome,split="\\.")[[1]][2])

bins <- binned.genome.reader(genome=genome, bin.size=20000, keep.rate=.9)
bins <- bins[[1]]
# remove unplaced chromosomes from bins
if(length(grep("_", bins$chr)) > 0){
  bins <- bins[-(grep("_", bins$chr)),]
}


bin.samp <- bins[sample(1:nrow(bins),10000),]

bin.sort = binSort(rep=rep, bins=bin.samp,TE.names=names(rep))
bin.samp = bin.sort$rates
bin.samp[,names(rep)] = bin.sort$rates[,names(rep)]/10000

repChoice = "L1PA"
rep.gr <- GRanges(seqnames = Rle(rep[[repChoice]]$genoName),
                  ranges = IRanges(start = rep[[repChoice]]$genoStart, end = rep[[repChoice]]$genoEnd))
bin.gr <- GRanges(seqnames=Rle(bin.samp$chr), 
                  ranges = IRanges(start=bin.samp$start, end=bin.samp$end-1))

diff.gr <- setdiff(bin.gr,rep.gr)
bin.diff <- as.matrix(findOverlaps(bin.gr, diff.gr))
if(length((1:nrow(bin.diff))[duplicated(bin.diff[,2])])!=0){
  print("some bins merged")
}
# bin diff provides our map for going back 



# we will get the GC percent of these regions and weight them by their length and get the mean

chromos <- unique(as.character(seqnames(diff.gr)))
bin.sortGC <- rep(NA, length(diff.gr))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=start(diff.gr)[as.character(seqnames(diff.gr)) == chromos[i]], 
                       end=end(diff.gr)[as.character(seqnames(diff.gr)) == chromos[i]]
  )
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[as.character(seqnames(diff.gr)) == chromos[i]] <- CGcontent
  print(chromos[i])
}

getContent <- data.frame(binID = unique(bin.diff[,1]), GCcontent = NA)
for( i in 1:nrow(getContent)){
  w <- width(diff.gr)[bin.diff[getContent$binID[i] == bin.diff[,1],2]]
  getContent[i,"GCcontent"] <- weighted.mean(bin.sortGC[bin.diff[getContent$binID[i] == bin.diff[,1],2]], w=w)
}



### old way

chromos <- as.character(unique(bin.samp$chr))
bin.sortGC <- rep(NA, nrow(bin.samp))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=bin.samp$start[bin.samp$chr == chromos[i]], 
                       end=bin.samp$end[bin.samp$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[bin.samp$chr == chromos[i]] <- CGcontent
  print(chromos[i])
}



# this way we can look for changes





plot(bin.sortGC[bin.samp[,repChoice] > 0],bin.samp[,repChoice][bin.samp[,repChoice] > 0], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))



modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)






plot(getContent$GCcontent,bin.samp[,repChoice][getContent$binID], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))


modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)


# GC content is being interesting, I wonder if getting the content of the regions then taking the mean over all the bases the right approach

# we seem to be taking the avarge rate
# rather than the rate over a position. 

# also it might be worth calculating the non TE GC rate. 
# it would be quite easy # we just take the setdiff
# and use those coordinates and map it back to the bin



# we could look at the difference between the two approaches 
# diff in GC contnet 
# against GC content

plot(bin.samp[getContent$binID,repChoice],getContent$GCcontent / bin.sortGC[getContent$binID])

plot(bin.samp[getContent$binID,repChoice],getContent$GCcontent - bin.sortGC[getContent$binID])
plot(getContent$GCcontent,bin.sortGC[getContent$binID])


# maybe we need to think about a different function for capturing cg adjustments



# thinking about how to control some potential sources of confunding 
# There will always be the problem of the repeat sequence contributing to the GC content 


# it is probably best to find at what level is the system robust to this 

# the probability a repeat picked at random has surrounding gc content of X




###### OLD APPROACH


chromos <- as.character(unique(bin.samp$chr))
bin.sortGC <- rep(NA, nrow(bin.samp))

for(i in 1:length(chromos)){
  Seq.set=DNAStringSet(bsSpec[[chromos[i]]], 
                       start=bin.samp$start[bin.samp$chr == chromos[i]], 
                       end=bin.samp$end[bin.samp$chr == chromos[i]])
  bases=data.frame(alphabetFrequency(Seq.set,baseOnly=TRUE))
  CGcontent <- rowSums(bases[,c("G", "C")]) / rowSums(bases[,c("A", "C", "T", "G")])
  bin.sortGC[bin.samp$chr == chromos[i]] <- CGcontent
  print(chromos[i])
}


# maybe just sample 1000 or something 




plot(bin.sortGC[bin.samp[,repChoice] > 0],bin.samp[,repChoice][bin.samp[,repChoice] > 0], pch = 16, cex = .3, ylim = c(0,1), xlim = c(0,1))

modLM <- lm(bin.samp[,repChoice][bin.samp[,repChoice] > 0] ~ bin.sortGC[bin.samp[,repChoice] > 0])
abline(modLM, lty = 2, lwd = 2, col = 2)

modLO <- loess(bin.samp[,repChoice] ~ bin.sortGC)
lines(seq(0,1,.01), predict(modLO, newdata=seq(0,1,.01)), col = 2, lwd = 2)


predict(modLO, .8)


TEs <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intergenic_reps$counts,refgene=refgene,type="intergenic",genome=genome,repList=rep)


# so now we have the GC regression slope

plot(TEs_intergenic$prime3gc/TEs_intergenic$baseFreq3prime, type = "l", col = 2)
par(new=TRUE)
plot(TEs_intergenic$rawRepCov3/TEs_intergenic$baseFreq3prime, type = "l", yaxt ="n")

axis(4,at=seq(0,1,by=.01), labels=seq(0,1,by=.01))

# how do we get the expected TE level based on GC content 
plot(predict(modLO, TEs_intergenic$prime3gc/TEs_intergenic$baseFreq3prime), type = "l")

lines(TEs_intergenic$rawRepCov3/TEs_intergenic$baseFreq3prime, col=2)




#TEs_intorn <- covCalcPlot5prime3primeGC(lenChoice=100000,repChoice="L1PA",repBins=intron_reps$counts,refgene=refgene,type="intron",genome=genome,repList=rep)



hist(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime)

plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, type = "l", col = 2)
par(new=TRUE)
plot(TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, type = "l", yaxt ="n")

axis(4,at=seq(0,1,by=.01), labels=seq(0,1,by=.01))

# how do we get the expected TE level based on GC content 
plot((1:100001),predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime), type = "l", ylim = c(0, .1))

lines((1:100001),TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, col=2)


# why do the values never go above a certain number



plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, TEs$rawRepCov3/TEs$baseFreq3prime)

plot(TEs$rawRepCov3/TEs$baseFreq3prime , predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime))



tail(predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime))

# why is it the expected value does not go above a certain limit
# is it becuase the GC range is out of the scope of the regression function

plot(TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime, type = "l")


gc.no <- TEs_intorn$prime3gc
# lets start doing the rolling mean appraoch of the gc content about 10000bp

# maybe it should be a rolling weighted mean
# one easy way would be to set up the coordinates 
# fill in the weights at each stage 
# fill in the values
# at each point we have the GC number or rate and the base frequency

library(zoo)

CG_coorect <- rollsum()


# we have a weighted mean problem 



z <- zoo(c(NA, 2, NA, 1, 4, 5, 2, NA))
na.fill(z)
na.fill(z, "extend")
na.fill(z, c("extend", NA))
na.fill(z, -(1:3))
na.fill(z, list(NA, NULL, NA))


x.Date <- as.Date(paste(2004, rep(1:4, 4:1), sample(1:28, 10), sep = "-"))
x <- zoo(rnorm(12), x.Date)

x
rollsum(x, 3,align="left" ,fill="extend")
rollsum(x, 3,align="right" ,fill="extend")

k = 20000
lenChoice = 100000


names(gc.no) <- 1:(lenChoice+1)

k = 1000
GCsum <- rollsum(TEs_intorn$prime3gc,k = k, fill="extend", align = "left")
baseSum <- rollsum(TEs_intorn$baseFreq3prime,k = k, fill="extend", align = "left")
TEsum <- rollsum(TEs_intorn$rawRepCov3,k = k, fill="extend", align = "left")
gcRate <- GCsum/baseSum


plot(TEsum/baseSum, type = "l")
lines(predict(modLO, gcRate), col = 2)

p = TEs_intorn$totalRep/TEs_intorn$totalCov
n = TEs_intorn$baseFreq3prime


# how would we adjust by GC content 

#lines((p*n)/n - predict(modLO, gcRate) + ((p*n)/n))
plot((TEsum/baseSum) - predict(modLO, gcRate) + ((p*n)/n) , type = "l")
lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)   , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  , col = 2  )

### we need to sample and heatmap this thing to see the dynamics



lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )

par(new=TRUE)
plot(n, type = "l")



lines( (p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n   )


##### the next thing is to normalise aginst all bins 
#### that is the best way to adhust the regression line 
#### furthermore we know this will change the CG dynamics



# why are we still getting these increases right near the front.

# we want p + the resdidual 




k = 1000


GCsum <- rollsum(TEs$prime3gc,k = k, fill="extend", align = "left")
baseSum <- rollsum(TEs$baseFreq3prime,k = k, fill="extend", align = "left")
TEsum <- rollsum(TEs$rawRepCov3,k = k, fill="extend", align = "left")
gcRate <- GCsum/baseSum


plot(TEsum/baseSum, type = "l", ylim=c(0, max(TEsum/baseSum)))
lines(predict(modLO, gcRate), col = 2)
lines(n*p/n)


p = TEs$totalRep/TEs$totalCov
n = TEs$baseFreq3prime


# how would we adjust by GC content 

#lines((p*n)/n - predict(modLO, gcRate) + ((p*n)/n))
plot((TEs$rawRepCov3/TEs$baseFreq3prime) - predict(modLO, gcRate) + ((p*n)/n) , type = "l")
lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)   , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n), col = 2  )





lines( ((p*n + 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )
lines( ((p*n - 2*((p*n) / sqrt( p*n*(1-p) )))/n)  - predict(modLO, gcRate) + ((p*n)/n) , col = 2  )



# work on correct sampling stratergies. 









plot(TEs_intorn$rawRepCov3/TEs_intorn$baseFreq3prime, type = "l")
lines(predict(modLO, gcRate),type = "l", col = 2)
lines(predict(modLO, TEs_intorn$prime3gc/TEs_intorn$baseFreq3prime), col = 3)

coords <- data.frame(start = seq(0-k, lenChoice, by = 1), end = seq(2, lenChoice + 1 +k +1, by = 1))
coords$start[coords$start < 1] <- 1
coords$end[coords$end > lenChoice + 1] <- lenChoice +1




A <- gc.no[coords$start : coords$end]




plot(coords$start,type = "l")
lines(coords$end)



