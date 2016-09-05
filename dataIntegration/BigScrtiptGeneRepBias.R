# lets put it all together 
# claculating the association with 






setwd("~/Desktop/Domain_manuscript/")

rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")
load(file = "~/Desktop/Domain_manuscript/R_objects/chromStateCombined")

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


# do an overlap to work out how many 

#121313687

intergenicJoinedFam <- data.frame(intergenic_reps$counts[,c("chr", "start", "end","Known")],
                                  old_L1 = rowSums(intergenic_reps$counts[,c("L1ME", "L1MD", "L1MC", "L1MB")]),
                                  new_L1 = rowSums(intergenic_reps$counts[,c("L1MA", "L1PB",  "L1PA","L1HS")]),
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
                new_L1 = rbind(rep$L1MA, rep$L1PB  ,rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))





intergenicSizesMean <- NULL
for(l in seq(0,max(intergenicJoinedFam$Known)+100,100)){
  intergenicSizesMean <- c(intergenicSizesMean, 
                           10^(mean(log10((intergenicJoinedFam$end - intergenicJoinedFam$start)[(intergenicJoinedFam$end - intergenicJoinedFam$start) > l]))))
}


intronSizesMean <- NULL
for(l in seq(0,max(intronJoinedFam$Known)+100,100)){
  intronSizesMean <- c(intronSizesMean, 10^(mean(log10((intronJoinedFam$end - intronJoinedFam$start)[(intronJoinedFam$end - intronJoinedFam$start) > l]))))
}

regions = c("intergenic", "intron")
TEcols <- c("red", "purple", "darkgreen","darkblue")



for(i in 1:4){
  for(j in 1:2){
    pdf(paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/biasLine/", regions[j],"/",names(joinRep)[i],".pdf", sep = "" ), onefile = T, height = 5, width = 5)
    
    TEfam = names(joinRep)[i]
    TEjoinFam = get(paste(regions[j], "JoinedFam", sep = ""))
    
    portion = TEjoinFam[,TEfam]/TEjoinFam$Known
    cuts <- cut(log10(TEjoinFam$Known), breaks = seq(2,8,by = .05))
    agg <- aggregate(x = portion,by = list(cuts), FUN = sum)
    # divide each number by how many times you see that cut
    summ <- table(cuts)
    
    n = sum(TEjoinFam[,TEfam])
    aggp <- aggregate(x = TEjoinFam$Known,by = list(cuts), FUN = sum)
    p = aggp$x/sum(aggp$x)
    m <- n*p
    SD <- sqrt(n*p*(1-p))
    
    plot(seq(2,8,by = .05)[summ > 0],agg$x/summ[summ > 0], main = "", xlim = c(1.9,7), 
         xlab = "interval length (kb)", ylab = "retrotransposon density", ylim = c(0,.22), pch = 16, col = TEcols[i], xaxt = "n")
    axis(side = 1, at = c(seq(0,8,by = 1)), labels = c((10^seq(0,8,by = 1))/1000))
    lines(seq(2,8,by = .05)[summ > 0],m/aggp$x, lty = 2)
    lines(seq(2,8,by = .05)[summ > 0],(m + (3*SD))/aggp$x)
    lines(seq(2,8,by = .05)[summ > 0],(m - (3*SD))/aggp$x)
    lo <- loess(agg$x/summ[summ > 0] ~ seq(2,8,by = .05)[summ > 0], span = .2)
    pred <- predict(lo,newdata = seq(2,8,by = .05),se = T)
    lines(seq(2,8,by = .05), pred$fit, col = TEcols[i], lwd = 3)
    # now we know what the preferance looks like we just have to apply it to real data. 
    yPoly = c(pred$fit + (3*pred$se.fit), rev(pred$fit - (3*pred$se.fit)))
    xPoly = c(seq(2,8,by = .05), seq(8,2,by = -.05))
    polygon(xPoly[complete.cases(yPoly)],yPoly[complete.cases(yPoly)], density = 20, col = TEcols[i], border = 0)
    
    assign(x = paste(regions[j], TEfam,  "Lo", sep = ""), value = lo)
    
    sizesMean <- get(paste(regions[j], "SizesMean", sep = ""))
    
    assign(x = paste(regions[j], TEfam,  "Pred", sep = ""), value = data.frame(position = (seq(0,max(TEjoinFam$Known)+100,100)/2),
                                                                     proportion = (predict(lo,newdata = log10(sizesMean))))
           )
    
    assign(x = paste(regions[j], TEfam,  "Bias", sep = ""), value = data.frame(position = (seq(0,max(TEjoinFam$Known)+100,100)/2),
                                                                     bias = (predict(lo,newdata = log10(sizesMean)))/(m/aggp$x)[1])
           )
    
    dev.off()
  }
}



intergenicLenChoice = max(intergenicJoinedFam$end - intergenicJoinedFam$start)/2
intronLenChoice = max(intronJoinedFam$end - intronJoinedFam$start)/2


region = "intron"
#TEfam = "new_L1"


for(i in 1:4){
  TEfam = names(joinRep)[i]
  lenChoice = get(paste(region, "LenChoice", sep =  ""))
  posStatBins = get(paste(region,"JoinedFam",sep = ""))
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")

  rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
  bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
  bpFreq[is.na(bpFreq)] <- 1
  
  rate <- rawCov/bpFreq
  

  
  # we can get both stats, 
  cutSite <- unique(as.integer(10^(seq(0,as.integer(log10(lenChoice+1))+1,.05))))
  cuted <- cut(1:(lenChoice+1), breaks = cutSite,right = F,ordered_result = T)
  cutPos<- as.integer(apply(data.frame(cutSite[1:(length(cutSite)-1)], cutSite[2:(length(cutSite))] ), 1, mean))
  Tab <- table(cuted)
  
  
  
  aggTEmean <- aggregate(rate, list(as.integer(cuted)), mean, simplify = F)
  aggTEsd <- aggregate(rate, list(as.integer(cuted)), sd, simplify = F)
  
  aggBPfreq <- aggregate(bpFreq, list(as.integer(cuted)), sum, simplify = F)
  aggBPraw <- aggregate(rawCov, list(as.integer(cuted)), sum, simplify = F)
  
  usePos <- cutPos[aggTEsd$Group.1]

  biasData <- get(paste(region, TEfam, "Bias", sep =""))
  SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
  biasPred <- predict(SPbiasData, usePos)
  
  shape.x <- c(log10(usePos), rev(log10(usePos)))
  shape.y <- c(aggTEmean$x + (2* aggTEsd$x), rev(aggTEmean$x - (2* aggTEsd$x)))
  shape.y.adj <- c((aggTEmean$x + (2* aggTEsd$x))/biasPred$y, rev((aggTEmean$x - (2* aggTEsd$x))/biasPred$y))
  
  sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))
  
  pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/corrected/", region, "/", TEfam,".pdf", sep = ""), height = 3, width = 10)
  par(mar = c(5,5,5,2))
  plot(log10(1:(lenChoice +1)), rate, col = 3, type = "n", ylim = c(0,.25), 
       main = "", ylab = "retrotransposon\ndensity", 
       xlab = "distance from boundary (log10 bp)")
  lines(log10(usePos), aggTEmean$x, type = "l", col = "grey60", lwd = 3)
  polygon(shape.x, shape.y, density =40,border = 0,col = "grey60")

  lines(log10(usePos), (aggBPraw$x/aggBPfreq$x)/biasPred$y, col = TEcols[i],lwd = 3 )
  polygon(shape.x, shape.y.adj, density =40,border = 0,col = TEcols[i],angle = 135)
  
  lines(log10(1:(lenChoice+1)), TEs_posStats$p + 2*sd, col = 1, lwd = 2)
  lines(log10(1:(lenChoice+1)), TEs_posStats$p - 2*sd, col = 1, lwd = 2)  
  legend("topleft",legend = c("uncorrected", "corrected"), fill = c("grey60", TEcols[i]), bty = "n")
  # question still is, are we accuratly capturing the bias effect through our smoothing
  dev.off()
  
  assign(x = paste(TEfam,region,"adjValues", sep = "_"),value = data.frame(position = usePos, rate = aggTEmean$x, adjRate = aggTEmean$x/biasPred$y)) 
  
}


# should we read in the data here ?




domainRanges <- read.table("Data/ConsTimingDomains", header = T)
domainRangesE <- domainRanges[domainRanges$domain == "ERD",]
domainRangesL <- domainRanges[domainRanges$domain == "LRD",]


domainE.gr <- GRanges(seqnames = Rle(domainRangesE$chr), 
                     ranges = IRanges(start = domainRangesE$start, end = domainRangesE$end))
domainE.gr

domainL.gr <- GRanges(seqnames = Rle(domainRangesL$chr), 
                     ranges = IRanges(start = domainRangesL$start, end = domainRangesL$end))
domainL.gr

findOverlaps(domainE.gr,domainL.gr)


# so now our domains are red into the data 

Intergenic.GR <- GRanges(seqnames = Rle(intergenic_reps$counts$chr),
                         ranges = IRanges(start = intergenic_reps$counts$start,
                                          end = intergenic_reps$counts$end)
)

OLintergenicE <- as.matrix(findOverlaps(Intergenic.GR, domainE.gr, type = "within"))
IntergenicE <- Intergenic.GR[OLintergenicE[,1]]

OLintergenicL <- as.matrix(findOverlaps(Intergenic.GR, domainL.gr, type = "within"))
IntergenicL <- Intergenic.GR[OLintergenicL[,1]]



plot(density(log10(width(IntergenicL))), ylim = c(0,.6), col = 2, xlab = "interval size (log10 bp)", main = "Intergenic regions in different genomic compartments")
lines(density(log10(width(IntergenicE))), col = 1)
legend("topright", legend = c("Early", "Late"), fill = c(1,2), title = "Replication Domain")

 
histL <- hist(log10(width(IntergenicL)), breaks = 100)
loL <- loess(histL$counts ~ histL$mids)
histE <- hist(log10(width(IntergenicE)), breaks = histL$breaks)
loE <- loess(histE$counts ~ histE$mids)
plot(seq(2,7,.01), predict(loL, seq(2,7,.01)), type = "l", ylim = c(0,200))
lines(seq(2,7,.01), predict(loE, seq(2,7,.01)), col = 2)

plot(log10(sort(width(IntergenicE))), (length(IntergenicE):1), ylab = "number of ranges > range width",
     xlab = "range width (log10 bp)", type = "l", main = "Intergenic position frequency", xlim = c(2,7))

lines(log10(sort(width(IntergenicL))), (length(IntergenicL):1), col = 2)
legend("topright", legend = c("Early", "Late"), fill = c(1,2), title = "Replication Domain")



# we need to fix our calling of constituative regions 

# or will rolling mean do the same?


# how many base sin each postion within our domains 
# calculate base frequency 


plot((sort(width(IntergenicE))), (length(IntergenicE):1), ylab = "number of ranges > range width",
     xlab = "range width (log10 bp)", type = "l", main = "Intergenic position frequency")

for(i in c("adjRate", "rate")){
  tabRes <- data.frame(Alu = c(NA,NA), new_L1 = c(NA,NA))
  for(j in c("Alu", "new_L1")){
    
    
    rownames(tabRes) <- c("Early", "Late")
    
    adjValue <- get(paste(j,"intergenic_adjValues", sep = "_"))
    SPadj <- smooth.spline(adjValue$position, adjValue[,i], all.knots = T)
    
    SPlevelsE <- smooth.spline(x = (sort(width(IntergenicE)))/2, y = (length(IntergenicE):1)*2)
    SPpredLevelE <- predict(SPlevelsE,1:(max(width(IntergenicE))/2))
    SpadjPredE <- predict(SPadj, 1:(max(width(IntergenicE))/2))
    Esum <- sum(SPpredLevelE$y * SpadjPredE$y)/sum(width(IntergenicE))
    
    
    SPlevelsL <- smooth.spline(x = (sort(width(IntergenicL)))/2, y = (length(IntergenicL):1)*2)
    SPpredLevelL <- predict(SPlevelsL,1:(max(width(IntergenicL))/2))
    SpadjPredL <- predict(SPadj, 1:(max(width(IntergenicL))/2))
    Lsum <- sum(SPpredLevelL$y * SpadjPredL$y)/sum(width(IntergenicL))
    
    tabRes[,j] <- c(Esum,Lsum)
  }
  tabRes["delta",] =  abs(tabRes[1,] - tabRes[2,])
  assign(x = paste(i,"Res",sep = ""), value = tabRes)
}



weighted.mean(x = c(Esum, Lsum),w = c(sum(width(IntergenicE)), sum(width(IntergenicL))))

## looks like the differences in ranges has nothing to do with it

head(intergenic_reps$rates)

EintergenicCounts <- intergenicJoinedFam[OLintergenicE[,1],]
LintergenicCounts <- intergenicJoinedFam[OLintergenicL[,1],]

ObsRes <- data.frame( Alu = c(sum(EintergenicCounts$Alu/sum(EintergenicCounts$Known)),
                              sum(LintergenicCounts$Alu/sum(LintergenicCounts$Known))),
                      new_L1 = c(sum(EintergenicCounts$new_L1/sum(EintergenicCounts$Known)),
                                 sum(LintergenicCounts$new_L1/sum(LintergenicCounts$Known)))
                      )
rownames(ObsRes) <- c("Early", "Late")
ObsRes["delta", ] <- abs(ObsRes[1,] - ObsRes[2,])



sum(intergenicJoinedFam$Alu)/sum(intergenicJoinedFam$Known)
sum(intergenicJoinedFam$new_L1)/sum(intergenicJoinedFam$Known)


layout(mat = matrix(c(1), nrow = 1))
pdf(file = paste("plots/geneRep/ModelFitting/", region,"/explainingDifference.pdf", sep = ""))
barplot(height = c(rateRes["delta","new_L1"]/ObsRes["delta","new_L1"], 
                   adjRateRes["delta","new_L1"]/ObsRes["delta","new_L1"],
                   rateRes["delta","Alu"]/ObsRes["delta","Alu"], 
                   adjRateRes["delta","Alu"]/ObsRes["delta","Alu"]
                   ),
        main = "difference between early and late replicating domains", 
        ylab = "proportion of explained difference",
        names = c("new_L1\nuncorrected", "new_L1\ncorrected","Alu\nuncorrected","Alu\ncorrected"),
        col = c(TEcols[2], TEcols[2], TEcols[3], TEcols[3]))
dev.off()

pdf(file = paste("plots/geneRep/ModelFitting/", region,"/modleData.pdf", sep = ""), onefile = T,width = 3.5,height = 3.5)
plot.new()
grid.table(round(rateRes, digits = 3))
mtext(text = "uncorrected model values",side = 3)
plot.new()
grid.table(round(adjRateRes, digits = 3))
mtext(text = "corrected model values",side = 3)
plot.new()
grid.table(round(ObsRes, digits = 3))
mtext(text = "observed values",side = 3)
dev.off()





##########
##
#   This is the enhancer area
######


# lets get enhancers as a rep type object

enhancers <- read.table("Data/FANTOM5/permissive_enhancers.bed.txt")
colnames(enhancers) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
enhancer.gr <- GRanges(seqnames = Rle(enhancers$chr),
                       ranges = IRanges(start = enhancers$start, end = enhancers$end))
enhancerOL <- as.matrix(findOverlaps(Intergenic.GR, enhancer.gr))
intergenicEnhancer <- intergenicJoinedFam[,c("chr", "start", "end","Known")]
enhanAgg <- aggregate(x = width(enhancer.gr[enhancerOL[,2]]), by = list(enhancerOL[,1]), FUN = sum)
intergenicEnhancer[enhanAgg$Group.1,"enhancer"] <- enhanAgg$x
intergenicEnhancer$enhancer[is.na(intergenicEnhancer$enhancer)] <- 0
enhancerList = list(enhancer = data.frame(genoName = enhancers$chr, genoStart = enhancers$start, genoEnd = enhancers$end))


EnhanRate <- sum(intergenicEnhancer$enhancer, na.rm = T)/sum(intergenicEnhancer$Known)


EnhanStatBins = intergenicEnhancer
lenChoice = max(EnhanStatBins$end - EnhanStatBins$start)/2
Enhan_posStatALL = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "enhancer", 
                                           repBins = EnhanStatBins, repList = enhancerList, 
                                           refgene = refgene, type = "intergenic", repType = "repeats")


EnhanStatBins = intergenicEnhancer[OLintergenicE[,1],]
lenChoice = max(EnhanStatBins$end - EnhanStatBins$start)/2
Enhan_posStatE <- covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "enhancer", 
                                          repBins = EnhanStatBins, repList = enhancerList, 
                                          refgene = refgene, type = "intergenic", repType = "repeats")


EnhanStatBins = intergenicEnhancer[OLintergenicL[,1],]
lenChoice = max(EnhanStatBins$end - EnhanStatBins$start)/2
Enhan_posStatL <- covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "enhancer", 
                                          repBins = EnhanStatBins, repList = enhancerList, 
                                          refgene = refgene, type = "intergenic", repType = "repeats")




layout(matrix(c(1,2), nrow = 1))

len = length(Enhan_posStatALL$rawRepCov5)
plot((1:len)[len:1], (Enhan_posStatALL$rawRepCov5/Enhan_posStatALL$baseFreq5prime)[1:len]/EnhanRate, 
     type = "l", xlim = c(10000,0), ylim = c(0,6),
     main = paste("enhancer intergenic upstream"))

len = length(Enhan_posStatE$rawRepCov5)
lines((1:len)[len:1], (Enhan_posStatE$rawRepCov5/Enhan_posStatE$baseFreq5prime)[1:len]/EnhanRate, col = 2)

len = length(Enhan_posStatL$rawRepCov5)
lines((1:len)[len:1], (Enhan_posStatL$rawRepCov5/Enhan_posStatL$baseFreq5prime)[1:len]/EnhanRate, col = 3)

# its porbably important i find a cool way of smoothing over that takes acoutn of the 
# element bp and the bp frequency at each position. 
# get a binning stratergy out of that. 



plot((1:len)[len:1], (Enhan_posStatALL$rawRepCov5/sum(intergenicEnhancer$enhancer, na.rm = T))[1:len], 
     type = "l", xlim = c(10000,0),
     main = paste("enhancer intergenic upstream"))

len = length(Enhan_posStatE$rawRepCov5)
lines((1:len)[len:1], (Enhan_posStatE$rawRepCov5/sum(intergenicEnhancer$enhancer, na.rm = T))[1:len], col = 2)

len = length(Enhan_posStatL$rawRepCov5)
lines((1:len)[len:1], (Enhan_posStatL$rawRepCov5/sum(intergenicEnhancer$enhancer, na.rm = T))[1:len], col = 3)




### 3prime
len = length(Enhan_posStatALL$rawRepCov5)
plot((1:len), (Enhan_posStatALL$rawRepCov3/Enhan_posStatALL$baseFreq3prime), 
     type = "l", ylim = c(0,.03), xlim = c(0,3000),
     main = paste("enhancer intergenic upstream"))

len = length(Enhan_posStatE$rawRepCov3)
lines((1:len), (Enhan_posStatE$rawRepCov3/Enhan_posStatE$baseFreq3prime), col = 2)

len = length(Enhan_posStatL$rawRepCov3)
lines((1:len), (Enhan_posStatL$rawRepCov3/Enhan_posStatL$baseFreq3prime), col = 3)

# could it be something to do with the enhancer density surrounding genes in a particular region
# what we in th eL regions is that enhancer levels are low

# or if we could suggest L1 insertion is particularly damaging to enhancer structure. 





region = "intergenic"
#TEfam = "new_L1"


for(i in 1:4){
  
  
  
  TEfam = names(joinRep)[i]
  
  
  
  posStatBins = rbind(get(paste("E",region,"Counts",sep = "")), get(paste("L",region,"Counts",sep = "")))
  lenChoice = max(posStatBins$end - posStatBins$start)/2
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")
  
  rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
  bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
  bpFreq[is.na(bpFreq)] <- 1
  Allrate <- rawCov/bpFreq

  TEs_posStatsAll <- TEs_posStats
  

  posStatBins = get(paste("E",region,"Counts",sep = ""))
  lenChoice = max(posStatBins$end - posStatBins$start)/2
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")
  
  rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
  bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
  bpFreq[is.na(bpFreq)] <- 1
  Erate <- rawCov/bpFreq

  TEs_posStatsE <- TEs_posStats
  

  posStatBins = get(paste("L",region,"Counts",sep = ""))
  lenChoice = max(posStatBins$end - posStatBins$start)/2
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")
  
  rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
  bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
  bpFreq[is.na(bpFreq)] <- 1
  Lrate <- rawCov/bpFreq
  
  TEs_posStatsL <- TEs_posStats
  
#   layout(matrix(c(1,2), nrow = 1))
#   len <- length(TEs_posStatsAll$rawRepCov5)
#   plot(log10(1:len)[len:1], (TEs_posStatsAll$rawRepCov5/TEs_posStatsAll$baseFreq5prime)[1:len], type = "l", ylim = c(0,.15), xlim = c(log10(len),0),
#        main = paste(TEfam, "intergenic upstream"))
#   len <- length(TEs_posStatsE$rawRepCov5)
#   lines(log10(1:len)[len:1], (TEs_posStatsE$rawRepCov5/TEs_posStatsE$baseFreq5prime)[1:len], col = 2)
#   len <- length(TEs_posStatsL$rawRepCov5)
#   lines(log10(1:len)[len:1], (TEs_posStatsL$rawRepCov5/TEs_posStatsL$baseFreq5prime)[1:len], col = 3)
#   abline(h=.06)
#   abline(v=2.5)
#   
#   
#   len <- length(TEs_posStatsAll$rawRepCov5)
#   plot(log10(1:len), (TEs_posStatsAll$rawRepCov3/TEs_posStatsAll$baseFreq3prime), type = "l", ylim = c(0,.15), xlim = c(0,log10(len)),
#        main = paste(TEfam, "intergenic downstream"))
#   len <- length(TEs_posStatsE$rawRepCov3)
#   lines(log10(1:len), (TEs_posStatsE$rawRepCov3/TEs_posStatsE$baseFreq3prime), col = 2)
#   len <- length(TEs_posStatsL$rawRepCov3)
#   lines(log10(1:len), (TEs_posStatsL$rawRepCov3/TEs_posStatsL$baseFreq3prime), col = 3)
#   abline(h=.06)
#   abline(v=2.5)
#   
  
  
  xlim = 3000
  ylim = .3
  
  layout(matrix(c(1,2), nrow = 1))
  par(mar=c(5,5,5,7))
  len <- length(TEs_posStatsAll$rawRepCov5)
  plot((1:len)[len:1], (TEs_posStatsAll$rawRepCov5/TEs_posStatsAll$baseFreq5prime)[1:len], type = "l", ylim = c(0,ylim), xlim = c(xlim,0),
       xlab = "position", ylab = "TE coverage per bp", main = paste(TEfam, region, "upstream"))
  grid()
  len <- length(TEs_posStatsE$rawRepCov5)
  lines((1:len)[len:1], (TEs_posStatsE$rawRepCov5/TEs_posStatsE$baseFreq5prime)[1:len], col = 2)
  len <- length(TEs_posStatsL$rawRepCov5)
  lines((1:len)[len:1], (TEs_posStatsL$rawRepCov5/TEs_posStatsL$baseFreq5prime)[1:len], col = 3)
 # abline(h=.05)
  #abline(h=.13)
  #abline(v=7000)
 
  
  par(new = T)
  len = length(Enhan_posStatALL$rawRepCov5)
  plot((1:len)[len:1], (Enhan_posStatALL$rawRepCov5/Enhan_posStatALL$baseFreq5prime)[1:len]/EnhanRate,
       type = "l", ylim = c(0,6), xlim = c(xlim,0),lty = 3, xaxt = "n", yaxt = "n",
       xlab = "", ylab = "")
  axis(4, at= seq(0,6),labels = round(seq(0,6),digits = 3),lty = 3)
  mtext("fold enrichment", side = 4, line = 3)
  
  len = length(Enhan_posStatE$rawRepCov5)
  lines((1:len)[len:1], (Enhan_posStatE$rawRepCov5/Enhan_posStatE$baseFreq5prime)[1:len]/EnhanRate, col = 2, lty = 3)
  
  len = length(Enhan_posStatL$rawRepCov5)
  lines((1:len)[len:1], (Enhan_posStatL$rawRepCov5/Enhan_posStatL$baseFreq5prime)[1:len]/EnhanRate, col = 3, lty =3)
  
  
  
  len <- length(TEs_posStatsAll$rawRepCov5)
  plot((1:len), (TEs_posStatsAll$rawRepCov3/TEs_posStatsAll$baseFreq3prime), type = "l", ylim = c(0,ylim), xlim = c(0,xlim),
       xlab = "position", ylab = "TE coverage per bp", main = paste(TEfam, region, "downstream"))
  grid()
  len <- length(TEs_posStatsE$rawRepCov3)
  lines((1:len), (TEs_posStatsE$rawRepCov3/TEs_posStatsE$baseFreq3prime), col = 2)
  len <- length(TEs_posStatsL$rawRepCov3)
  lines((1:len), (TEs_posStatsL$rawRepCov3/TEs_posStatsL$baseFreq3prime), col = 3)
#   abline(h=.05)
#   abline(h=.13)
#   abline(v=7000)
  
  
  par(new = T)
  len = length(Enhan_posStatALL$rawRepCov3)
  plot((1:len), (Enhan_posStatALL$rawRepCov3/Enhan_posStatALL$baseFreq3prime)/EnhanRate,
       type = "l", ylim = c(0,6), xlim = c(0,xlim),lty = 3, xaxt = "n", yaxt = "n",
       ylab = "", xlab = "")
  axis(4, at= seq(0,6),labels = round(seq(0,6),digits = 3),lty = 3)
  mtext("enhancer fold enrichment", side = 4, line = 3)
  
  len = length(Enhan_posStatE$rawRepCov5)
  lines((1:len), (Enhan_posStatE$rawRepCov3/Enhan_posStatE$baseFreq3prime)/EnhanRate, col = 2, lty = 3)
  
  len = length(Enhan_posStatL$rawRepCov3)
  lines((1:len), (Enhan_posStatL$rawRepCov3/Enhan_posStatL$baseFreq3prime)/EnhanRate, col = 3, lty =3)
  
  
  
}



# this isnt making sense, for L1s we had the expectation that the sizes of regions
# were casuing a major differnce in the acculation of elements. 


# we say that bias is caused by the difference between open and closed chromatin mapping to different size locations. 



# I think we can see now at least there is likly very little boundary effect
# However there is strong regional effect. 



# something about the open state thats different than the closed

# so we need to figuer out why we are getting the enrichemnts of elements in differetn areas





files <- list.files("Data/OpenChromSynth/")
allOpen <- NULL
for(i in 1:length(files)){
  dat <- read.table(paste("Data/OpenChromSynth/",files[i], sep = ""))
  dat[,22] <- as.factor(gsub("_[0-9]+","",dat[,4]))
  colnames(dat) <- c("chr", "start", "end", "hitTypeNumber", "score", "strand", "start2", "end2", "colour", "Pval",
                     "dnaseSig", "dnasePval", "faireSig", "fairePval", "polIISig", "polIIPval", "ctcfSig","ctcfPval","cmycSig", "cmycPval", "ocCode", "ocType")
  name <- gsub(pattern = "_OC.bed", replacement = "",x = files[i])
  assign(name, dat)
  allOpen <- c(allOpen, list(get(name)))
}
names(allOpen) <- gsub(pattern = "_OC.bed", replacement = "", x = files)


allCell <- NULL
for(i in 2:length(allOpen)){
  OpenCell <- allOpen[[i]]
  OpenCell[,"cellName"] <- names(allOpen)[i]
  allCell <- rbind(allCell,OpenCell)
}



#OpenCell <- allOpen$H1_ES
#OpenCell <- OpenCell[OpenCell$ocType == "DNaseOnly",]


OpenCell <- allCell
OpenCell <- OpenCell[OpenCell$dnasePval > 6,]

Open.gr <- reduce(GRanges(seqnames = Rle(OpenCell$chr),
                       ranges = IRanges(start = OpenCell$start, end = OpenCell$end)))
Open.gr <- intersect(Open.gr,Intergenic.GR)
openOL <- as.matrix(findOverlaps(Intergenic.GR, Open.gr))
intergenicOpen <- intergenicJoinedFam[,c("chr", "start", "end","Known")]
openAgg <- aggregate(x = width(Open.gr[openOL[,2]]), by = list(openOL[,1]), FUN = sum)
intergenicOpen[openAgg$Group.1,"openChrom"] <- openAgg$x
intergenicOpen$openChrom[is.na(intergenicOpen$openChrom)] <- 0
OpenChromList = list(openChrom = data.frame(genoName = OpenCell$chr, genoStart = OpenCell$start, genoEnd = OpenCell$end))


OpenRate <- sum(intergenicOpen$openChrom, na.rm = T)/sum(intergenicOpen$Known)




OpenStatBins = intergenicOpen
# OpenStatBins$start=OpenStatBins$start-5000
# OpenStatBins$end = OpenStatBins$end + 5000
# OpenStatGR <- GRanges(seqnames = Rle(OpenStatBins$chr),
#                       ranges = IRanges(start = OpenStatBins$start, end = OpenStatBins$end))
#OpenStatBins = OpenStatBins[countOverlaps(OpenStatGR) == 1 & countOverlaps(query = OpenStatGR, subject = refgene.gr) == 2,]
lenChoice = max(OpenStatBins$end - OpenStatBins$start)/2
Open_posStatALL = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "openChrom", 
                                           repBins = OpenStatBins, repList = OpenChromList, 
                                           refgene = refgene, type = "intergenic", repType = "repeats")


OpenStatBins = intergenicOpen[OLintergenicE[,1],]
lenChoice = max(OpenStatBins$end - OpenStatBins$start)/2
Open_posStatE = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "openChrom", 
                                          repBins = OpenStatBins, repList = OpenChromList, 
                                          refgene = refgene, type = "intergenic", repType = "repeats")


OpenStatBins = intergenicOpen[OLintergenicL[,1],]
lenChoice = max(OpenStatBins$end - OpenStatBins$start)/2
Open_posStatL = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "openChrom", 
                                          repBins = OpenStatBins, repList = OpenChromList, 
                                          refgene = refgene, type = "intergenic", repType = "repeats")






i = 3
TEfam = names(joinRep)[i]


region = "intergenic"
posStatBins = rbind(get(paste("E",region,"Counts",sep = "")), get(paste("L",region,"Counts",sep = "")))
# posStatBins$start=posStatBins$start-5000
# posStatBins$end = posStatBins$end + 5000
# posStatGR <- GRanges(seqnames = Rle(posStatBins$chr),
#                       ranges = IRanges(start = posStatBins$start, end = posStatBins$end))
# posStatBins = posStatBins[countOverlaps(posStatGR) == 1 & countOverlaps(query = posStatGR, subject = refgene.gr) == 2,]
lenChoice = max(posStatBins$end - posStatBins$start)/2
TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")



rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
bpFreq[is.na(bpFreq)] <- 1
Allrate <- rawCov/bpFreq

TEs_posStatsAll <- TEs_posStats


posStatBins = get(paste("E",region,"Counts",sep = ""))
lenChoice = max(posStatBins$end - posStatBins$start)/2

TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")

rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
bpFreq[is.na(bpFreq)] <- 1
Erate <- rawCov/bpFreq

TEs_posStatsE <- TEs_posStats


posStatBins = get(paste("L",region,"Counts",sep = ""))
lenChoice = max(posStatBins$end - posStatBins$start)/2

TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")

rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
bpFreq <- (TEs_posStats$baseFreq3prime + TEs_posStats$baseFreq5prime[(lenChoice+1):1])
bpFreq[is.na(bpFreq)] <- 1
Lrate <- rawCov/bpFreq

TEs_posStatsL <- TEs_posStats






layout(matrix(c(1,2),nrow = 2))

lenChoice = 10000

plot((1:lenChoice), (TEs_posStatsAll$rawRepCov3/TEs_posStatsAll$baseFreq3prime)[1:lenChoice] , lty=1, , type="l",ylim = c(0,.3),)
lines((1:lenChoice), (TEs_posStatsE$rawRepCov3/TEs_posStatsE$baseFreq3prime)[1:lenChoice] , lty=1, col = 2)
lines((1:lenChoice), (TEs_posStatsL$rawRepCov3/TEs_posStatsL$baseFreq3prime)[1:lenChoice] , lty=1, col = 3)

par(new=T)
plot((1:lenChoice), (Open_posStatALL$rawRepCov3/Open_posStatALL$baseFreq3prime)[1:lenChoice] , lty = 3, type = "l", ylim = c(0,.3))
lines((1:lenChoice), (Open_posStatE$rawRepCov3/Open_posStatE$baseFreq3prime)[1:lenChoice] , col = 2, lty = 3)
lines((1:lenChoice), (Open_posStatL$rawRepCov3/Open_posStatL$baseFreq3prime)[1:lenChoice] , col = 3, lty = 3)

plot((1:lenChoice), rev(TEs_posStatsAll$rawRepCov5/TEs_posStatsAll$baseFreq5prime)[1:lenChoice] , lty=1, , type="l",ylim = c(0,.3))
lines((1:lenChoice), rev(TEs_posStatsE$rawRepCov5/TEs_posStatsE$baseFreq5prime)[1:lenChoice] , lty=1, col = 2)
lines((1:lenChoice), rev(TEs_posStatsL$rawRepCov5/TEs_posStatsL$baseFreq5prime)[1:lenChoice] , lty=1, col = 3)

par(new=T)
plot((1:lenChoice), rev(Open_posStatALL$rawRepCov5/Open_posStatALL$baseFreq5prime)[1:lenChoice] , lty = 3, type = "l", ylim = c(0,.3))
lines((1:lenChoice), rev(Open_posStatE$rawRepCov5/Open_posStatE$baseFreq5prime)[1:lenChoice] , col = 2, lty = 3)
lines((1:lenChoice),  rev(Open_posStatL$rawRepCov5/Open_posStatL$baseFreq5prime)[1:lenChoice] , col = 3, lty = 3)



ratioEL5Open <- rev(Open_posStatE$rawRepCov5/Open_posStatE$baseFreq5prime)[1:lenChoice]/rev(Open_posStatL$rawRepCov5/Open_posStatL$baseFreq5prime)[1:lenChoice]
ratioEL5Alu <- rev(TEs_posStatsE$rawRepCov5/TEs_posStatsE$baseFreq5prime)[1:lenChoice]/rev(TEs_posStatsL$rawRepCov5/TEs_posStatsL$baseFreq5prime)[1:lenChoice]

ratioEL3Open <- (Open_posStatE$rawRepCov3/Open_posStatE$baseFreq3prime)[1:lenChoice]/(Open_posStatL$rawRepCov3/Open_posStatL$baseFreq3prime)[1:lenChoice]
ratioEL3Alu <- (TEs_posStatsE$rawRepCov3/TEs_posStatsE$baseFreq3prime)[1:lenChoice]/(TEs_posStatsL$rawRepCov3/TEs_posStatsL$baseFreq3prime)[1:lenChoice]

plot(ratioEL5Alu, type = "l", ylim = c(0,5))
lines(ratioEL5Open, col= 1, lty = 3)

plot(ratioEL3Alu, type = "l", ylim = c(0,5))
lines(ratioEL3Open, col= 1, lty = 3)



# lets say the ratios don't match up. 
# Is this something to do with differences in the small regions

OpenStatBins = intergenicOpen[OLintergenicE[,1],]
plot(log10(OpenStatBins$Known), (OpenStatBins$openChrom/(OpenStatBins$Known)), pch = 16, cex = .3)

plot(log10(EintergenicCounts$Known),(EintergenicCounts$Alu/EintergenicCounts$Known), col = 2, cex = .1, pch = 16)
points(log10(LintergenicCounts$Known),(LintergenicCounts$Alu/LintergenicCounts$Known), col = 3, cex = .1, pch = 16)


OpenStatBins = intergenicOpen[OLintergenicE[,1],]
points(log10(OpenStatBins$Known), (OpenStatBins$openChrom/(OpenStatBins$Known)), pch = 16, cex = .1, col = 4)
OpenStatBins = intergenicOpen[OLintergenicL[,1],]
points(log10(OpenStatBins$Known), (OpenStatBins$openChrom/(OpenStatBins$Known)), pch = 16, cex = .1, col = 5)



layout(1)
OpenStatBins = intergenicOpen[OLintergenicL[,1],]

plot((OpenStatBins$openChrom)/LintergenicCounts$Known, (LintergenicCounts$new_L1)/LintergenicCounts$Known, pch = 16, cex = .3)

