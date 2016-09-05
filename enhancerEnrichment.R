# maybe we should look at somethign like enhacers and see what TE distribtuions are around them

# it would be interesting to see if it behaves a certain way




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

# from here we will get our good part of the genome our intergenic/intronic set where we can pull out clusters we want

nonExonRegion.gr <- GRanges(seqnames = Rle(c(as.character(bins_gene_gap$chr),as.character(bins_intron$chr))),
                            ranges = IRanges(start = c(bins_gene_gap$start,bins_intron$start),
                                             end = c(bins_gene_gap$end,bins_intron$end)))


web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chrom_info <- read.delim(textConnection(txt), header = FALSE)




rep = rep_info(spec1=spec1,genome=genome)


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB,rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))
insertSize = 1
TPRTinsertNew_L1 <- joinRep$new_L1
TPRTinsertNew_L1 <- TPRTinsertNew_L1[(TPRTinsertNew_L1$repStart < 0 & TPRTinsertNew_L1$repStart > -42) | (TPRTinsertNew_L1$repLeft < 0 & TPRTinsertNew_L1$repLeft > -42),]
TPRTinsertNew_L1$genoStart[TPRTinsertNew_L1$strand == "+"] <- TPRTinsertNew_L1$genoEnd[TPRTinsertNew_L1$strand == "+"] - insertSize
TPRTinsertNew_L1$genoEnd <- TPRTinsertNew_L1$genoStart + insertSize
joinRep$TPRT <- TPRTinsertNew_L1

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


# 
# 
# 
# enhancers <- read.table("Data/FANTOM5/permissive_enhancers.bed.txt")
# colnames(enhancers) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
# enhancer.gr <- GRanges(seqnames = Rle(enhancers$chr),
#                        ranges = IRanges(start = enhancers$start, end = enhancers$end))


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
for(i in 1:length(allOpen)){
  OpenCell <- allOpen[[i]]
  OpenCell[,"cellName"] <- names(allOpen)[i]
  allCell <- rbind(allCell,OpenCell)
}
openCells <- allCell[allCell$dnasePval > 1.3,]
openCellBumpy.gr <- GRanges(seqnames = Rle(openCells$chr),
                              ranges = IRanges(start = openCells$start, end = openCells$end))

openCellBumpy.gr <- subsetByOverlaps(openCellBumpy.gr, nonExonRegion.gr, type = "within")

openCell.gr <- reduce(openCellBumpy.gr)
#end(openCell.gr) <- start(openCell.gr)
#openCell.gr = openCell.gr[width(openCell.gr) > 1000 & width(openCell.gr) < 3000]
gapsOpenCellAll.gr <- gaps((openCell.gr))
gapsOpenCellAll.gr <- gapsOpenCellAll.gr[width(gapsOpenCellAll.gr) < 7000000]

#openCell.gr <- enhancer.gr

# gapsOpenCellAll.gr <- gaps((openCell.gr))
# gapsOpenCellAll.gr <- gapsOpenCellAll.gr[width(gapsOpenCellAll.gr) < 7000000]
# gapsOpenCellAll.gr <- gapsOpenCellAll.gr[NottooSmall]
gapsOpenCellE.gr <- subsetByOverlaps(gapsOpenCellAll.gr,domainE.gr, type = "within")
gapsOpenCellL.gr <- subsetByOverlaps(gapsOpenCellAll.gr,domainL.gr, type = "within")



# some of the base stats for open cells 


openCellsE.gr <- subsetByOverlaps(openCell.gr, domainE.gr, type = "within")
openCellsL.gr <- subsetByOverlaps(openCell.gr, domainL.gr, type = "within")
openCellsAll.gr <- openCell.gr


layout(1)

Eactive <- hist(countOverlaps(openCellsE.gr,openCellBumpy.gr), breaks = seq(0,200), freq = F)
Lactive <- hist(countOverlaps(openCellsL.gr,openCellBumpy.gr), breaks = seq(0,200), freq = F)

EnonExonBases.gr <- intersect(nonExonRegion.gr,domainE.gr) 
LnonExonBases.gr <- intersect(nonExonRegion.gr,domainL.gr) 

pdf(file = "plots/TEopenChromInteract/clusters.pdf", height = 3.5,width = 2.5)
plot(Eactive$mids + .5, log2((Eactive$counts/sum(width(EnonExonBases.gr)))* (1*10^7)), col = 2,xlim = c(1,50), type = "l", 
     xlab = "DNase1 HS per cluster", ylab = "DNase1 clusters per 10 Mb", lwd = 3, yaxt = "n", ylim = c(-10,10))
axis(side = 2, at = seq(-10,10, 5), labels = round(2^(seq(-10,10, 5)),2))
lineG <- log2((Lactive$counts/sum(width(LnonExonBases.gr)))* (1*10^7))
lineG[is.infinite(lineG)] <- -10
lines(Lactive$mids + .5,lineG, col = 3, lwd = 3)
legend("topright", legend = c("cERD", "cLRD"),title = "Domain", fill = c(2,3), bty = "n", cex = .7)
dev.off()
# this shows how active each region is 

# next is to look at what is happeing inside the gaps and what is happening outside

# lets just compare the mean to the total 
AllcovTE <- LcovTE <- EcovTE <-matrix(data = NA, nrow = 4,ncol = 2)

for(i in 1:4){
  lenChoice = c(2000,10000,3000,3000, 5000)[i]
  repChoice = c("Alu", "new_L1", "old_L1", "Ancient")[i]
  
  TE.gr <- GRanges(seqnames = Rle(joinRep[[repChoice]]$genoName), 
                   ranges = IRanges(start = joinRep[[repChoice]]$genoStart, end = joinRep[[repChoice]]$genoEnd))
  
  nonExonTE.gr <- intersect(TE.gr,nonExonRegion.gr)
  clusterTE.gr <- intersect(TE.gr, openCellsAll.gr)
  AllcovTE[i,] = rev(c(sum(width(clusterTE.gr))/sum(width(openCellsAll.gr)), 
                   sum(width(nonExonTE.gr))/sum(as.numeric(width(nonExonRegion.gr)))))
  nonExonTE.gr <- intersect(TE.gr,EnonExonBases.gr)
  clusterTE.gr <- intersect(TE.gr, openCellsE.gr)
  EcovTE[i,] = rev(c(sum(width(clusterTE.gr))/sum(width(openCellsE.gr)), 
                   sum(width(nonExonTE.gr))/sum(as.numeric(width(EnonExonBases.gr)))))
  nonExonTE.gr <- intersect(TE.gr,LnonExonBases.gr)
  clusterTE.gr <- intersect(TE.gr, openCellsL.gr)
  LcovTE[i,] = rev(c(sum(width(clusterTE.gr))/sum(width(openCellsL.gr)), 
                   sum(width(nonExonTE.gr))/sum(as.numeric(width(LnonExonBases.gr)))))

}
rownames(EcovTE) <- rownames(LcovTE) <- c("Alu", "new_L1", "old_L1", "Ancient")

pdf(file = "plots/TEopenChromInteract/TEOLcluster.pdf", height = 7,width = 3.3)
layout(c(1,2))
barplot(t(EcovTE),beside = T,space = c(0,.3), ylim = c(0,.23), 
        col = c(rep("darkgreen",2), rep("purple", 2), rep("red", 2), rep("darkblue", 2)), 
        density = c(-1,30), xaxt = "n", las = 2)
legend("topright", legend = c("non exon", "DNase cluster"), density= c(-1,30), box.lwd = 0)
barplot(t(LcovTE), beside = T, ylim = c(0,.23),space = c(0,.3), 
        col = c(rep("darkgreen",2), rep("purple", 2), rep("red", 2), rep("darkblue", 2)), 
        density = c(-1,30), las = 2)
dev.off()


# uncorrected levels
layout(matrix(1:4,nrow = 2))

layout(1)
i = 5
#for(i in 1:4){
  
  lenChoice = c(2000,10000,3000,3000, 5000)[i]
  repChoice = c("Alu", "new_L1", "old_L1", "Ancient", "TPRT")[i]
  

  
  TE.gr <- GRanges(seqnames = Rle(joinRep[[repChoice]]$genoName), 
                   ranges = IRanges(start = joinRep[[repChoice]]$genoStart, end = joinRep[[repChoice]]$genoEnd))
  
  gapsOpenCellAll <- data.frame(chr = seqnames(gapsOpenCellAll.gr), start = start(gapsOpenCellAll.gr), end = end(gapsOpenCellAll.gr), Known = width(gapsOpenCellAll.gr))
  gapsOpenCellAll[,repChoice] = countOverlaps(gapsOpenCellAll.gr, TE.gr)
  resAll <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellAll,repList = joinRep)
  #plot(c(resAll$rawRepCov5/resAll$baseFreq5,resAll$rawRepCov3/resAll$baseFreq3), type = "l")
  
  gapsOpenCellE<- data.frame(chr = seqnames(gapsOpenCellE.gr), start = start(gapsOpenCellE.gr), end = end(gapsOpenCellE.gr), Known = width(gapsOpenCellE.gr))
  gapsOpenCellE[,repChoice] = countOverlaps(gapsOpenCellE.gr, TE.gr)
  resE <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellE,repList = joinRep)
  #plot(c(res$rawRepCov5/res$baseFreq5,res$rawRepCov3/res$baseFreq3), ylim = c(0,.2), type = "l")
  
  gapsOpenCellL<- data.frame(chr = seqnames(gapsOpenCellL.gr), start = start(gapsOpenCellL.gr), end = end(gapsOpenCellL.gr), Known = width(gapsOpenCellL.gr))
  gapsOpenCellL[,repChoice] = countOverlaps(gapsOpenCellL.gr, TE.gr)
  resL <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellL,repList = joinRep)
  #plot(c(res$rawRepCov5/res$baseFreq5,res$rawRepCov3/res$baseFreq3), ylim = c(0,.2), type = "l")
  
  
  
  plot(c(resAll$rawRepCov5/resAll$baseFreq5,resAll$rawRepCov3/resAll$baseFreq3), ylim = c(0,.001), type = "h", xaxt = "n",
       xlab = "distance from enhancer (bp)", ylab = "repeat coverage", main = repChoice)
  axis(side = 1,at = c(seq(lenChoice, 0,length.out = 5),seq(lenChoice, lenChoice*2,length.out = 5)), 
       labels = c(seq(0, lenChoice,length.out = 5),seq(0, lenChoice,length.out = 5)))
  lines(c(resE$rawRepCov5/resE$baseFreq5,resE$rawRepCov3/resE$baseFreq3), col = 2)
  lines(c(resL$rawRepCov5/resL$baseFreq5,resL$rawRepCov3/resL$baseFreq3), col = 3)
  
#}


# choose TE here

AllSizesMean <- NULL
ESizesMean <- NULL
LSizesMean <- NULL
for(l in seq(0,20000,10)){
  AllSizesMean <- c(AllSizesMean, 
                           10^(mean(log10(width(gapsOpenCellAll.gr)[width(gapsOpenCellAll.gr) > l]))))
  ESizesMean <- c(ESizesMean, 
                  10^(mean(log10(width(gapsOpenCellE.gr)[width(gapsOpenCellE.gr) > l]))))
  LSizesMean <- c(LSizesMean, 
                  10^(mean(log10(width(gapsOpenCellL.gr)[width(gapsOpenCellL.gr) > l]))))
}



for(i in 1:5){
  TEcols = c("darkgreen", "purple", "red", "darkblue", "orange")[i]
  repChoice = c("Alu", "new_L1", "old_L1", "Ancient", "TPRT")[i]
  TE.gr <- GRanges(seqnames = Rle(joinRep[[repChoice]]$genoName), 
                   ranges = IRanges(start = joinRep[[repChoice]]$genoStart, end = joinRep[[repChoice]]$genoEnd))
  # choose gap region here 
  for(j in 1:3){
    
    gType = c("All","E", "L")[j]
    
    pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/biasLine/dnase1/", repChoice,gType,"OpendnaseNonExon.pdf", sep = ""), height = 5, width = 5)

    gapDat.gr <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
    gapTable <- data.frame(chr = seqnames(gapDat.gr), start = start(gapDat.gr), end = end(gapDat.gr), Known = width(gapDat.gr))
    
    TE.int <- intersect(TE.gr,gapDat.gr)
    OL <- as.matrix(findOverlaps(TE.int, gapDat.gr))
    agg <- aggregate(width(TE.int)[OL[,1]], by = list(OL[,2]), sum)
    
    gapTable[agg$Group.1, "TE"] <- agg$x
    gapTable$TE[is.na(gapTable$TE)] <- 0
    
    
    portion = gapTable$TE/gapTable$Known
    
    
    cuts <- cut(log10(gapTable$Known), breaks = seq(1,8,by = .05))
    agg <- aggregate(x = portion,by = list(cuts), FUN = sum)
    agg <- merge(x = data.frame(Group.1 = levels(cuts)), y = agg, by.x = 1, by.y = 1,all=TRUE)
    agg$x[is.na(agg$x)] <- 0
    # divide each number by how many times you see that cut
    summ <- table(cuts)
    
    n = sum(gapTable$TE)
    aggp <- aggregate(x = gapTable$Known,by = list(cuts), FUN = sum)
    aggp <- merge(x = data.frame(Group.1 = levels(cuts)), y = aggp, by.x = 1, by.y = 1,all=TRUE)
    aggp$x[is.na(aggp$x)] <- 0
    
    
    p = aggp$x/sum(aggp$x)
    m <- n*p
    SD <- sqrt(n*p*(1-p))
    
    seqNo <- seq(1,7.95,by = .05)[summ > 0]
    plot(seqNo,agg$x[summ > 0]/summ[summ > 0], main = "", xlim = c(.9,7), 
         xlab = "interval length (kb)", ylab = "retrotransposon density", ylim = c(0,.3), pch = 16, col = TEcols, xaxt = "n")
    axis(side = 1, at = c(seq(0,8,by = 1)), labels = c((10^seq(0,8,by = 1))/1000))
    lines(seqNo,m[summ > 0]/aggp$x[summ > 0], lty = 2)
    lines(seqNo,(m + (3*SD))[summ > 0]/aggp$x[summ > 0])
    lines(seqNo,(m - (3*SD))[summ > 0]/aggp$x[summ > 0])
    lo <- loess(agg$x[summ > 0]/summ[summ > 0] ~ seqNo,span = .2)
#    pred <- predict(lo,newdata = seq(2,8,by = .05),se = T)
    pred <- predict(lo,newdata = seqNo, se = T)
#    lines(seq(2,8,by = .05), pred$fit, col = TEcols[i], lwd = 3)
    lines(seqNo, pred$fit, col = TEcols, lwd = 3)
    
    yPoly = c(pred$fit + (3*pred$se.fit), rev(pred$fit - (3*pred$se.fit)))
#    xPoly = c(seq(2,8,by = .05), seq(8,2,by = -.05))
    xPoly = c(seqNo, rev(seqNo))
    polygon(xPoly[complete.cases(yPoly)],yPoly[complete.cases(yPoly)], density = 20, col = TEcols, border = 0)
    
    
    assign(x = paste(gType, repChoice,  "Lo", sep = ""), value = lo)
    
    sizesMean <- get(paste(gType, "SizesMean", sep = ""))
    
    assign(x = paste(gType, repChoice,  "Pred", sep = ""), value = data.frame(position = (seq(0,20000,10)/2),
                                                                               proportion = (predict(lo,newdata = log10(sizesMean))))
    )
    
    assign(x = paste(gType, repChoice,  "Bias", sep = ""), value = data.frame(position = (seq(0,20000,10)/2),
                                                                               bias = (predict(lo,newdata = log10(sizesMean)))/(m[summ > 0]/aggp$x[summ > 0])[1])
    )
    
    dev.off()
  }
  
}








# so we can establish that there is reduced coverage in enhancer regions
# would be interesting to see enhancer density in each intronic region. 
# it is most interesting to see what is happening with different sizes


### maybe we could extend this to look at TFBSs and other regulatory sites.  



TEcols = c("darkgreen", "purple", "red", "darkblue", "orange")
for(j in 1:3){
for(i in 1:4){
  TEfam = c("Alu", "new_L1", "old_L1", "Ancient", "TPRT")[i]
  lenChoice = c(2000,10000,3000,3000, 5000)[i]
  gType = c("All", "E", "L")[j]
  
  TE.gr <- GRanges(seqnames = Rle(joinRep[[TEfam]]$genoName), 
                   ranges = IRanges(start = joinRep[[TEfam]]$genoStart, end = joinRep[[TEfam]]$genoEnd))
  
  gapsOpen.gr  <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
  
  gapsOpen <- data.frame(chr = seqnames(gapsOpen.gr), start = start(gapsOpen.gr), end = end(gapsOpen.gr), Known = width(gapsOpen.gr))
  gapsOpen[,TEfam] = countOverlaps(gapsOpen.gr, TE.gr)
  
  
  posStatBins = gapsOpen
  
  TEs_posStats <- covCalcPlot(lenChoice = lenChoice, repChoice = TEfam, repBins = posStatBins,repList = joinRep)
  
  rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
  bpFreq <- (TEs_posStats$baseFreq3 + TEs_posStats$baseFreq5[(lenChoice+1):1])
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
  
  biasData <- get(paste(gType, TEfam, "Bias", sep =""))
  SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
  biasPred <- predict(SPbiasData, usePos)
  

  shape.x <- c((usePos), rev((usePos)))
  shape.y <- c(aggTEmean$x + (2* aggTEsd$x), rev(aggTEmean$x - (2* aggTEsd$x)))
  shape.y.adj <- c((aggTEmean$x + (2* aggTEsd$x))/biasPred$y, rev((aggTEmean$x - (2* aggTEsd$x))/biasPred$y))
  
  sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))
  
  pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/corrected/dnase1/", gType, TEfam,".pdf", sep = ""), height = 3, width = 10)
  par(mar = c(5,5,5,2))
  plot((1:(lenChoice +1)), rate, col = 3, type = "n", ylim = c(0,.35), 
       main = "", ylab = "retrotransposon\ndensity", 
       xlab = "distance from boundary (bp)", yaxt = "n")
  axis(side = 2,at = seq(0,1,.1), las = 2)
  lines((usePos), aggTEmean$x, type = "l", col = "grey60", lwd = 3)
  polygon(shape.x, shape.y, density =40,border = 0,col = "grey60")
  
  lines((usePos), (aggBPraw$x/aggBPfreq$x)/biasPred$y, col = TEcols[i],lwd = 3 )
  polygon(shape.x, shape.y.adj, density =40,border = 0,col = TEcols[i],angle = 135)
  
  lines((1:(lenChoice+1)), TEs_posStats$p + 2*sd, col = 1, lwd = 2)
  lines((1:(lenChoice+1)), TEs_posStats$p - 2*sd, col = 1, lwd = 2)  
  legend("topright",legend = c("uncorrected", "corrected"), fill = c("grey60", TEcols[i]), bty = "n")
  # question still is, are we accuratly capturing the bias effect through our smoothing
  dev.off()
  
  #assign(x = paste(TEfam,region,"adjValues", sep = "_"),value = data.frame(position = usePos, rate = aggTEmean$x, adjRate = aggTEmean$x/biasPred$y)) 
  
}

}



# we can rerun this loop to get the lines we need
pdf(file = "plots/TEopenChromInteract/TEboundary.pdf", height = 5, width =3)
layout(c(1,2,3,4))
par(mar = c(2,5,.5,5))
for(i in 1:4){
   TEfam = c("Alu", "new_L1", "old_L1", "Ancient")[i]
    lenChoice = 12000
    xlimChoice = c(2000,8000,2000,2000)[i]
  plot((1:(lenChoice +1)), (1:(lenChoice +1)), col = 3, type = "n", ylim = c(-.1,.3), xlim = c(0, xlimChoice), 
         ylab = "", xlab = "", xaxt = "n", las = 1, yaxt = "n")
  mtext(TEfam,side = 4,line = 1)
  axis(side = 1,at = c(0,(xlimChoice/2),xlimChoice),labels = c(0,(xlimChoice/2)/1000,xlimChoice/1000))
  axis(side = 2,at = c(-.1, 0, .1, .2,.3), las = 2)
  
  for(j in 1:2){
   
    gType = c("E", "L")[j]
    
    TE.gr <- GRanges(seqnames = Rle(joinRep[[TEfam]]$genoName), 
                     ranges = IRanges(start = joinRep[[TEfam]]$genoStart, end = joinRep[[TEfam]]$genoEnd))
    
    gapsOpen.gr  <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
    
    gapsOpen <- data.frame(chr = seqnames(gapsOpen.gr), start = start(gapsOpen.gr), end = end(gapsOpen.gr), Known = width(gapsOpen.gr))
    gapsOpen[,TEfam] = countOverlaps(gapsOpen.gr, TE.gr)
    
    posStatBins = gapsOpen
    
    TEs_posStats <- covCalcPlot(lenChoice = lenChoice, repChoice = TEfam, repBins = posStatBins,repList = joinRep)
    
    rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
    bpFreq <- (TEs_posStats$baseFreq3 + TEs_posStats$baseFreq5[(lenChoice+1):1])
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
    
    biasData <- get(paste(gType, TEfam, "Bias", sep =""))
    SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
    biasPred <- predict(SPbiasData, usePos)
    
    sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))

    shape.x <- c((usePos), rev((usePos)))
    shape.y <- c((((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p) + (3 *sd[usePos]), 
                 rev((((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p) - (3 * sd[usePos])))

    polygon(shape.x, shape.y, density =20,border = 0,col = j+1 ,angle = 135)

    lines((usePos), (((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p), col = j + 1,lwd = 1 )

  abline(h = 0, lty = 2)   

  }
  
}

dev.off()




##### the TRPT stuff


lenChoice = c(30000)
repChoice = c("TPRT")



TE.gr <- GRanges(seqnames = Rle(joinRep[[repChoice]]$genoName), 
                 ranges = IRanges(start = joinRep[[repChoice]]$genoStart, end = joinRep[[repChoice]]$genoEnd))

gapsOpenCellAll <- data.frame(chr = seqnames(gapsOpenCellAll.gr), start = start(gapsOpenCellAll.gr), end = end(gapsOpenCellAll.gr), Known = width(gapsOpenCellAll.gr))
gapsOpenCellAll[,repChoice] = countOverlaps(gapsOpenCellAll.gr, TE.gr)
resAll <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellAll,repList = joinRep)
#plot(c(resAll$rawRepCov5/resAll$baseFreq5,resAll$rawRepCov3/resAll$baseFreq3), type = "l")

gapsOpenCellE<- data.frame(chr = seqnames(gapsOpenCellE.gr), start = start(gapsOpenCellE.gr), end = end(gapsOpenCellE.gr), Known = width(gapsOpenCellE.gr))
gapsOpenCellE[,repChoice] = countOverlaps(gapsOpenCellE.gr, TE.gr)
resE <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellE,repList = joinRep)
#plot(c(res$rawRepCov5/res$baseFreq5,res$rawRepCov3/res$baseFreq3), ylim = c(0,.2), type = "l")

gapsOpenCellL<- data.frame(chr = seqnames(gapsOpenCellL.gr), start = start(gapsOpenCellL.gr), end = end(gapsOpenCellL.gr), Known = width(gapsOpenCellL.gr))
gapsOpenCellL[,repChoice] = countOverlaps(gapsOpenCellL.gr, TE.gr)
resL <- covCalcPlot(lenChoice = lenChoice, repChoice = repChoice, repBins = gapsOpenCellL,repList = joinRep)
#plot(c(res$rawRepCov5/res$baseFreq5,res$rawRepCov3/res$baseFreq3), ylim = c(0,.2), type = "l")




pdf(file = "plots/TEopenChromInteract/TPRT.pdf", height = 4, width = 3)

layout(c(1))
cutter <- cut(sqrt(1:(lenChoice+1)), breaks = seq(0, sqrt(lenChoice),by = 10 ),include.lowest = T)
raw3 = aggregate(x = resE$rawRepCov3, list(cutter), FUN = sum)
freq3 = aggregate(x = resE$baseFreq3, list(cutter), FUN = sum)
raw5 = aggregate(x = rev(resE$rawRepCov5), list(cutter), FUN = sum)
freq5 = aggregate(x = rev(resE$baseFreq5), list(cutter), FUN = sum)

plot(seq(5,by = 10,length.out = nrow(raw5) )^2,((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000, 
     ylim = c(0,350), type = "l", col = 2, lwd = 3, xlab = "distance from DNase clusters",
     ylab = "new_L1 insertion per Mb")
abline(h = (sum(width(intersect(TE.gr,EnonExonBases.gr)))/sum(width(EnonExonBases.gr)))* 10^6 , lty = 2, lwd = 3, col=2)
#axis(side = 1,at = seq(5,by = 10,length.out = nrow(raw5)*10 )^2)

cutter <- cut(sqrt(1:(lenChoice+1)), breaks = seq(0, sqrt((lenChoice+1)),by = 10 ),include.lowest = T)
raw3 = aggregate(x = resL$rawRepCov3, list(cutter), FUN = sum)
freq3 = aggregate(x = resL$baseFreq3, list(cutter), FUN = sum)
raw5 = aggregate(x = rev(resL$rawRepCov5), list(cutter), FUN = sum)
freq5 = aggregate(x = rev(resL$baseFreq5), list(cutter), FUN = sum)

lines(seq(5,by = 10,length.out = nrow(raw5) )^2,((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000, 
      ylim = c(0,350), type = "l",lwd = 3, col = 3)
abline(h = (sum(width(intersect(TE.gr,LnonExonBases.gr)))/sum(width(LnonExonBases.gr)))* 10^6 , lty = 2, lwd = 3, col=3)

dev.off()





# pick a range of criteria, so FAIRE and DNase

# If we could generate some insertion sites for L1s, if we use 3' end mapping to look at where they stack up

# merge our data from a range of cell types



# some of the background stuff, TE density per intergenic region



# produce gap data here 


# differences between types of dnase hypersenstivity sites

