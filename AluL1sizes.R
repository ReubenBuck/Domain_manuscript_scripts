#### L1/alu size analysis analysis 

# we need three figs for this part 
# the coverage levels 

# the position ratio


# O-E 


setwd("~/Desktop/Domain_manuscript/")

rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")


rep = rep_info(spec1=spec1,genome=genome)


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB  ,rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2))




domainRanges <- read.table("Data/ConsTimingDomains", header = T)
domainRangesE <- domainRanges[domainRanges$domain == "ERD",]
domainRangesL <- domainRanges[domainRanges$domain == "LRD",]

domainE.gr <- GRanges(seqnames = Rle(domainRangesE$chr), 
                      ranges = IRanges(start = domainRangesE$start, end = domainRangesE$end))
domainE.gr

domainL.gr <- GRanges(seqnames = Rle(domainRangesL$chr), 
                      ranges = IRanges(start = domainRangesL$start, end = domainRangesL$end))
domainL.gr





new_L1repstart <- joinRep$new_L1$repStart
new_L1repstart[joinRep$new_L1$strand == "-"] <- joinRep$new_L1$repLeft[joinRep$new_L1$strand == "-"]
new_L1repstart[joinRep$new_L1$repEnd - new_L1repstart < 0] <- joinRep$new_L1$repEnd[joinRep$new_L1$repEnd - new_L1repstart < 0]

new_L1repleft <- joinRep$new_L1$repLeft
new_L1repleft[joinRep$new_L1$strand == "-"] <- joinRep$new_L1$repStart[joinRep$new_L1$strand == "-"]
new_L1repleft <- abs(new_L1repleft)

L1Range <- GRanges(seqnames = joinRep$new_L1$genoName, 
                     ranges = IRanges(start = joinRep$new_L1$genoStart, end = joinRep$new_L1$genoEnd),
                     repRange = IRanges(end = (new_L1repleft + joinRep$new_L1$repEnd) - new_L1repstart, start = (new_L1repleft + joinRep$new_L1$repEnd) - joinRep$new_L1$repEnd))


#### here we can begin to compare
Alurepstart <- joinRep$Alu$repStart
Alurepstart[joinRep$Alu$strand == "-"] <- joinRep$Alu$repLeft[joinRep$Alu$strand == "-"]
Alurepstart[joinRep$Alu$repEnd - Alurepstart < 0] <- joinRep$Alu$repEnd[joinRep$Alu$repEnd - Alurepstart < 0]

Alurepleft <- joinRep$Alu$repLeft
Alurepleft[joinRep$Alu$strand == "-"] <- joinRep$Alu$repStart[joinRep$Alu$strand == "-"]
Alurepleft <- abs(Alurepleft)

AluRange <- GRanges(seqnames = joinRep$Alu$genoName, 
                   ranges = IRanges(start = joinRep$Alu$genoStart, end = joinRep$Alu$genoEnd),
                   repRange = IRanges(end = (Alurepleft + joinRep$Alu$repEnd) - Alurepstart, start = (Alurepleft + joinRep$Alu$repEnd) - joinRep$Alu$repEnd))



AluE <- subsetByOverlaps(AluRange, domainE.gr)
AluL <- subsetByOverlaps(AluRange, domainL.gr)
AluEcov <- coverage(elementMetadata(AluE)$repRange)[1:300]
AluLcov <- coverage(elementMetadata(AluL)$repRange)[1:300]
AluAcov <- coverage(elementMetadata(AluRange)$repRange)[1:300]


L1E <- subsetByOverlaps(L1Range, domainE.gr)
L1L <- subsetByOverlaps(L1Range, domainL.gr)
L1Ecov <- coverage(elementMetadata(L1E)$repRange)[1:6000]
L1Lcov <- coverage(elementMetadata(L1L)$repRange)[1:6000]
L1Acov <- coverage(elementMetadata(L1Range)$repRange)[1:6000]


pdf(file = "~/Desktop/Domain_manuscript/plots/TEopenChromInteract/RTNpos.pdf", height = 3.2,width = 3)
par(mar = c(5,5,2,2))
plot((1:6000),(L1Acov/sum(L1Acov))*6000,type = "l", ylim = c(0,4), col = "purple", lwd = 3, xaxt = "n", xlab = "", 
     ylab = "", las = 2)
par(new = TRUE)
plot((AluAcov/sum(AluAcov))*300, type = "l", ylim = c(0,4), col = "darkgreen", lwd = 3,axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(0,150,300), labels = FALSE)
mtext(text = c(0,50,100),side = 1,line = 1,at = c(0,150,300), col = "black")
mtext(text = "position in full length RTN\nfrom 3' end (%)",side = 1,line = 3)
mtext(text = "RTN position density (O/E)", side = 2,line = 2)
abline(h = 1, lty = 3, lwd = 3)
legend("topright", legend = c("Alu", "new L1"), fill = c("darkgreen", "purple"), bty = "n")
dev.off()
# mtext(text = c(1,3000,6000),side = 1,line = 1,at = c(1,150,300), col = "purple")
# mtext(text = c(1,150,300),side = 1,line = 2,at = c(1,150,300), col = "darkgreen")

# to pos 300
plot((L1Acov/sum(L1Acov))*6000,type = "l", ylim = c(0,4), xlim=c(0,300),col = "purple", lwd = 3, xaxt = "n", xlab = "", ylab = "")
lines((AluAcov/sum(AluAcov))*300, type = "l", col = "darkgreen", lwd = 3)
axis(side = 1, at = c(1,150,300))




pdf(file = "~/Desktop/Domain_manuscript/plots/TEopenChromInteract/RTNdomainEnrich.pdf", height = 3.2,width = 3)
par(mar = c(5,5,2,2))
plot((L1Ecov/sum(width(L1E))) / (L1Lcov/sum(width(L1L))),type = "l", ylim = c(.8,1.5), col = "purple", lwd = 3, xaxt = "n", 
     xlab = "", ylab = "", las = 2)
par(new = TRUE)
plot((AluEcov/sum(width(AluE))) / (AluLcov/sum(width(AluL))), type = "l", 
     ylim = c(.8,1.5), col = "darkgreen", lwd = 3,axes = FALSE, xlab = "", ylab = "")
axis(side = 1, at = c(1,150,300), labels = FALSE)
mtext(text = c(0,50,100),side = 1,line = 1,at = c(0,150,300), col = "black")
mtext("position in full length RTN\nfrom 3' end (%)",side = 1,line = 3)
mtext("RTN position density ratio\n(cERD:cLRD)", side = 2, line = 3)
abline(h = 1, lty = 3, lwd = 3)
dev.off()

# mtext(text = c(1,3000,6000),side = 1,line = 1,at = c(1,150,300), col = "purple")
# mtext(text = c(1,150,300),side = 1,line = 2,at = c(1,150,300), col = "darkgreen")


### the next bit is to measure the corrected TPRT O-E cERD:cLRD ratio


###### read in all our region defining stuff

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
gapsOpenCellAll.gr <- gaps((openCell.gr))
gapsOpenCellAll.gr <- gapsOpenCellAll.gr[width(gapsOpenCellAll.gr) < 7000000]


gapsOpenCellE.gr <- subsetByOverlaps(gapsOpenCellAll.gr,domainE.gr, type = "within")
gapsOpenCellL.gr <- subsetByOverlaps(gapsOpenCellAll.gr,domainL.gr, type = "within")





insertSize = 1
TPRTinsertNew_L1 <- joinRep$new_L1
TPRTinsertNew_L1 <- TPRTinsertNew_L1[(TPRTinsertNew_L1$repStart < 0 & TPRTinsertNew_L1$repStart > -50) | (TPRTinsertNew_L1$repLeft < 0 & TPRTinsertNew_L1$repLeft > -50),]
TPRTinsertNew_L1$genoStart[TPRTinsertNew_L1$strand == "+"] <- TPRTinsertNew_L1$genoEnd[TPRTinsertNew_L1$strand == "+"] - insertSize
TPRTinsertNew_L1$genoEnd <- TPRTinsertNew_L1$genoStart + insertSize
joinRep$new_L1_TPRT <- TPRTinsertNew_L1

TPRTinsertAlu <- joinRep$Alu
TPRTinsertAlu <- TPRTinsertAlu[(TPRTinsertAlu$repStart < 0 & TPRTinsertAlu$repStart > -50) | (TPRTinsertAlu$repLeft < 0 & TPRTinsertAlu$repLeft > -50),]
TPRTinsertAlu$genoStart[TPRTinsertAlu$strand == "+"] <- TPRTinsertAlu$genoEnd[TPRTinsertAlu$strand == "+"] - insertSize
TPRTinsertAlu$genoEnd <- TPRTinsertAlu$genoStart + insertSize
joinRep$Alu_TPRT <- TPRTinsertAlu





EnonExonBases.gr <- intersect(nonExonRegion.gr,domainE.gr) 
LnonExonBases.gr <- intersect(nonExonRegion.gr,domainL.gr) 


for(i in 1:2){
  
  lenChoice = c(30000)
  repChoice = c("Alu_TPRT", "new_L1_TPRT")[i]
  
  
  
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
  
  
  
  #pdf(file = "plots/TEopenChromInteract/TPRT.pdf", height = 4, width = 3)
  
  # if we could simply correct the insertion rate then it would be easy to do what we need
  
  
  layout(c(1))
  cutter <- cut(sqrt(1:(lenChoice+1)), breaks = seq(0, sqrt(lenChoice),by = 10 ),include.lowest = T)
  raw3 = aggregate(x = resE$rawRepCov3, list(cutter), FUN = sum)
  freq3 = aggregate(x = resE$baseFreq3, list(cutter), FUN = sum)
  raw5 = aggregate(x = rev(resE$rawRepCov5), list(cutter), FUN = sum)
  freq5 = aggregate(x = rev(resE$baseFreq5), list(cutter), FUN = sum)
  
  plot(seq(5,by = 10,length.out = nrow(raw5) )^2,((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000, 
       ylim = c(0,2000), type = "l", col = 2, lwd = 3, xlab = "distance from DNase clusters",
       ylab = "new_L1 insertion per Mb")
  abline(h = (sum(width(intersect(TE.gr,EnonExonBases.gr)))/sum(width(EnonExonBases.gr)))* 10^6 , lty = 2, lwd = 3, col=2)
  #axis(side = 1,at = seq(5,by = 10,length.out = nrow(raw5)*10 )^2)
  ob <- (((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000)
  ex <- ((sum(width(intersect(TE.gr,EnonExonBases.gr)))/sum(width(EnonExonBases.gr)))* 10^6)
  
  assign(x = paste("E",repChoice, "OE", sep = ""),
         value = ob/ex)
  
  
  
  cutter <- cut(sqrt(1:(lenChoice+1)), breaks = seq(0, sqrt((lenChoice+1)),by = 10 ),include.lowest = T)
  raw3 = aggregate(x = resL$rawRepCov3, list(cutter), FUN = sum)
  freq3 = aggregate(x = resL$baseFreq3, list(cutter), FUN = sum)
  raw5 = aggregate(x = rev(resL$rawRepCov5), list(cutter), FUN = sum)
  freq5 = aggregate(x = rev(resL$baseFreq5), list(cutter), FUN = sum)
  
  lines(seq(5,by = 10,length.out = nrow(raw5) )^2,((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000, 
        ylim = c(0,350), type = "l",lwd = 3, col = 3)
  abline(h = (sum(width(intersect(TE.gr,LnonExonBases.gr)))/sum(width(LnonExonBases.gr)))* 10^6 , lty = 2, lwd = 3, col=3)
  
  ob <- ((raw3$x+raw5$x) / (freq3$x + freq5$x )) * 1000000
  ex <- (sum(width(intersect(TE.gr,LnonExonBases.gr)))/sum(width(LnonExonBases.gr)))* 10^6
  assign(x = paste("L",repChoice, "OE", sep = ""),
         value = ob/ex)
  
}

# so lets work on this when we get back

#dev.off()



# I don't think this trick works because it starts to water down the effects
# maybe look at values relative to an expected out come related to a lowerbound sigma
pdf(file = "~/Desktop/Domain_manuscript/plots/TEopenChromInteract/RTNdomainDNase1.pdf", height = 3.2,width = 3)
par(mar = c(5,5,2,2))
plot(seq(5,by = 10,length.out = nrow(raw5) )^2,(EAlu_TPRTOE) / (LAlu_TPRTOE), type = "l",ylim = c(0,8), xlim = c(0,15000),col = "darkgreen", lwd = 3, 
     ylab = "", xlab = "distance from DNase1\ncluster boundary (kb)", xaxt = "n", las = 2)
lines(seq(5,by = 10,length.out = nrow(raw5) )^2,(Enew_L1_TPRTOE)/(Lnew_L1_TPRTOE), lwd = 3, col = "purple")
abline(h=1, lty = 3, lwd = 3)
mtext("O/E RTN insertion ratio\n(cERD:cLRD)",side = 2,line = 2)
axis(side = 1,at = seq(0,30000,5000), seq(0,30,5))
dev.off()

plot(seq(5,by = 10,length.out = nrow(raw5) )^2,EAlu_TPRTOE, type = "l", ylim = c(0,10))
lines(seq(5,by = 10,length.out = nrow(raw5) )^2,LAlu_TPRTOE, lty = 2)

lines(seq(5,by = 10,length.out = nrow(raw5) )^2,Enew_L1_TPRTOE, col = 2)
lines(seq(5,by = 10,length.out = nrow(raw5) )^2,Lnew_L1_TPRTOE, col = 2, lty = 3)

abline(h=1)

# maybe there's some sort of transformation I can apply



#### lets try this in the morning


# 
# 
# 
# 
# AllSizesMean <- NULL
# ESizesMean <- NULL
# LSizesMean <- NULL
# for(l in seq(0,20000,10)){
#   AllSizesMean <- c(AllSizesMean, 
#                     10^(mean(log10(width(gapsOpenCellAll.gr)[width(gapsOpenCellAll.gr) > l]))))
#   ESizesMean <- c(ESizesMean, 
#                   10^(mean(log10(width(gapsOpenCellE.gr)[width(gapsOpenCellE.gr) > l]))))
#   LSizesMean <- c(LSizesMean, 
#                   10^(mean(log10(width(gapsOpenCellL.gr)[width(gapsOpenCellL.gr) > l]))))
# }
# 
# 
# 
# for(i in 1:2){
#   TEcols = c("darkgreen", "purple")[i]
#   repChoice = c("Alu_TPRT", "new_L1_TPRT")[i]
#   TE.gr <- GRanges(seqnames = Rle(joinRep[[repChoice]]$genoName), 
#                    ranges = IRanges(start = joinRep[[repChoice]]$genoStart, end = joinRep[[repChoice]]$genoEnd))
#   # choose gap region here 
#   for(j in 1:3){
#     
#     gType = c("All","E", "L")[j]
#     
#     pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/biasLine/dnase1/", repChoice,gType,"OpendnaseNonExon.pdf", sep = ""), height = 5, width = 5)
#     
#     gapDat.gr <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
#     gapTable <- data.frame(chr = seqnames(gapDat.gr), start = start(gapDat.gr), end = end(gapDat.gr), Known = width(gapDat.gr))
#     
#     TE.int <- intersect(TE.gr,gapDat.gr)
#     OL <- as.matrix(findOverlaps(TE.int, gapDat.gr))
#     agg <- aggregate(width(TE.int)[OL[,1]], by = list(OL[,2]), sum)
#     
#     gapTable[agg$Group.1, "TE"] <- agg$x
#     gapTable$TE[is.na(gapTable$TE)] <- 0
#     
#     
#     portion = gapTable$TE/gapTable$Known
#     
#     cuts <- cut(sqrt(gapTable$Known), breaks = seq(0, 2000,by = 10 ),include.lowest = T)
#     #cuts <- cut(log10(gapTable$Known), breaks = seq(1,8,by = .1))
#     agg <- aggregate(x = portion,by = list(cuts), FUN = sum)
#     agg <- merge(x = data.frame(Group.1 = levels(cuts)), y = agg, by.x = 1, by.y = 1,all=TRUE)
#     agg$x[is.na(agg$x)] <- 0
#     # divide each number by how many times you see that cut
#     summ <- table(cuts)
#     
#     n = sum(gapTable$TE)
#     aggp <- aggregate(x = gapTable$Known,by = list(cuts), FUN = sum)
#     aggp <- merge(x = data.frame(Group.1 = levels(cuts)), y = aggp, by.x = 1, by.y = 1,all=TRUE)
#     aggp$x[is.na(aggp$x)] <- 0
#     
#     
#     p = aggp$x/sum(aggp$x)
#     m <- n*p
#     SD <- sqrt(n*p*(1-p))
#     
#     seqNo <- seq(0, 2000 - 10,by = 10 )[summ > 0]
#     plot(seqNo,agg$x[summ > 0]/summ[summ > 0], main = "", xlim = c(0,2000), 
#          xlab = "interval length (kb)", ylab = "RTN density", ylim = c(0,.002), pch = 16, col = TEcols, xaxt = "n")
#     axis(side = 1, at = c(seq(0,2000,by = 100)), labels = seq(0,2000,by = 100)^2)
#     lines(seqNo,m[summ > 0]/aggp$x[summ > 0], lty = 2)
#     lines(seqNo,(m + (3*SD))[summ > 0]/aggp$x[summ > 0])
#     lines(seqNo,(m - (3*SD))[summ > 0]/aggp$x[summ > 0])
#     lo <- loess(agg$x[summ > 0]/summ[summ > 0] ~ seqNo,span = .2)
#     #    pred <- predict(lo,newdata = seq(2,8,by = .05),se = T)
#     pred <- predict(lo,newdata = seqNo, se = T)
#     #    lines(seq(2,8,by = .05), pred$fit, col = TEcols[i], lwd = 3)
#     lines(seqNo, pred$fit, col = TEcols, lwd = 3)
#     
#     yPoly = c(pred$fit + (3*pred$se.fit), rev(pred$fit - (3*pred$se.fit)))
#     #    xPoly = c(seq(2,8,by = .05), seq(8,2,by = -.05))
#     xPoly = c(seqNo, rev(seqNo))
#     polygon(xPoly[complete.cases(yPoly)],yPoly[complete.cases(yPoly)], density = 20, col = TEcols, border = 0)
#     
#     
#     assign(x = paste(gType, repChoice,  "Lo", sep = ""), value = lo)
#     
#     sizesMean <- get(paste(gType, "SizesMean", sep = ""))
#     
#     assign(x = paste(gType, repChoice,  "Pred", sep = ""), value = data.frame(position = (seq(0,20000,10)/2),
#                                                                               proportion = (predict(lo,newdata = log10(sizesMean))))
#     )
#     
#     assign(x = paste(gType, repChoice,  "Bias", sep = ""), value = data.frame(position = (seq(0,20000,10)/2),
#                                                                               bias = (predict(lo,newdata = log10(sizesMean)))/(m[summ > 0]/aggp$x[summ > 0])[1])
#     )
#     
#     dev.off()
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # so we can establish that there is reduced coverage in enhancer regions
# # would be interesting to see enhancer density in each intronic region. 
# # it is most interesting to see what is happening with different sizes
# 
# 
# ### maybe we could extend this to look at TFBSs and other regulatory sites.  
# 
# 
# 
# TEcols = c("darkgreen", "purple")
# for(j in 1:3){
#   for(i in 1:2){
#     TEfam = c("Alu_TPRT", "new_L1_TPRT")[i]
#     lenChoice = c(2000,2000)[i]
#     gType = c("All", "E", "L")[j]
#     
#     TE.gr <- GRanges(seqnames = Rle(joinRep[[TEfam]]$genoName), 
#                      ranges = IRanges(start = joinRep[[TEfam]]$genoStart, end = joinRep[[TEfam]]$genoEnd))
#     
#     gapsOpen.gr  <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
#     
#     gapsOpen <- data.frame(chr = seqnames(gapsOpen.gr), start = start(gapsOpen.gr), end = end(gapsOpen.gr), Known = width(gapsOpen.gr))
#     gapsOpen[,TEfam] = countOverlaps(gapsOpen.gr, TE.gr)
#     
#     
#     posStatBins = gapsOpen
#     
#     TEs_posStats <- covCalcPlot(lenChoice = lenChoice, repChoice = TEfam, repBins = posStatBins,repList = joinRep)
#     
#     rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
#     bpFreq <- (TEs_posStats$baseFreq3 + TEs_posStats$baseFreq5[(lenChoice+1):1])
#     bpFreq[is.na(bpFreq)] <- 1
#     
#     rate <- rawCov/bpFreq
#     
#     # we can get both stats, 
#     cutSite <- unique(as.integer(10^(seq(0,as.integer(log10(lenChoice+1))+1,.05))))
#     cuted <- cut(1:(lenChoice+1), breaks = cutSite,right = F,ordered_result = T)
#     cutPos<- as.integer(apply(data.frame(cutSite[1:(length(cutSite)-1)], cutSite[2:(length(cutSite))] ), 1, mean))
#     Tab <- table(cuted)
#     
#     aggTEmean <- aggregate(rate, list(as.integer(cuted)), mean, simplify = F)
#     aggTEsd <- aggregate(rate, list(as.integer(cuted)), sd, simplify = F)
#     
#     aggBPfreq <- aggregate(bpFreq, list(as.integer(cuted)), sum, simplify = F)
#     aggBPraw <- aggregate(rawCov, list(as.integer(cuted)), sum, simplify = F)
#     
#     usePos <- cutPos[aggTEsd$Group.1]
#     
#     biasData <- get(paste(gType, TEfam, "Bias", sep =""))
#     SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
#     biasPred <- predict(SPbiasData, usePos)
#     
#     
#     shape.x <- c((usePos), rev((usePos)))
#     shape.y <- c(aggTEmean$x + (2* aggTEsd$x), rev(aggTEmean$x - (2* aggTEsd$x)))
#     shape.y.adj <- c((aggTEmean$x + (2* aggTEsd$x))/biasPred$y, rev((aggTEmean$x - (2* aggTEsd$x))/biasPred$y))
#     
#     sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))
#     
#     pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/corrected/dnase1/", gType, TEfam,".pdf", sep = ""), height = 3, width = 10)
#     plot((1:(lenChoice +1)), rate, col = 3, type = "n", ylim = c(0,.002), 
#          main = "", ylab = "RTN density", 
#          xlab = "distance from boundary (bp)", yaxt = "n")
#     axis(side = 2,at = seq(0,1,.1), las = 2)
#     lines((usePos), aggTEmean$x, type = "l", col = "grey60", lwd = 3)
#     polygon(shape.x, shape.y, density =40,border = 0,col = "grey60")
#     
#     lines((usePos), (aggBPraw$x/aggBPfreq$x)/biasPred$y, col = TEcols[i],lwd = 3 )
#     polygon(shape.x, shape.y.adj, density =40,border = 0,col = TEcols[i],angle = 135)
#     
#     lines((1:(lenChoice+1)), TEs_posStats$p + 2*sd, col = 1, lwd = 2)
#     lines((1:(lenChoice+1)), TEs_posStats$p - 2*sd, col = 1, lwd = 2)  
#     legend("topright",legend = c("uncorrected", "corrected"), fill = c("grey60", TEcols[i]), bty = "n")
#     # question still is, are we accuratly capturing the bias effect through our smoothing
#     dev.off()
#     
#     #assign(x = paste(TEfam,region,"adjValues", sep = "_"),value = data.frame(position = usePos, rate = aggTEmean$x, adjRate = aggTEmean$x/biasPred$y)) 
#     
#   }
#   
# }
# 
# 
# 
# # we can rerun this loop to get the lines we need
# pdf(file = "plots/TEopenChromInteract/TPRT_TEboundary.pdf", height = 5, width =3)
# layout(c(1,2,3,4))
# par(mar = c(2,5,.5,5))
# for(i in 1:4){
#   TEfam = c("Alu", "new_L1", "old_L1", "Ancient")[i]
#   lenChoice = 12000
#   xlimChoice = c(2000,8000,2000,2000)[i]
#   plot((1:(lenChoice +1)), (1:(lenChoice +1)), col = 3, type = "n", ylim = c(-.1,.3), xlim = c(0, xlimChoice), 
#        ylab = "", xlab = "", xaxt = "n", las = 1, yaxt = "n")
#   mtext(TEfam,side = 4,line = 1)
#   axis(side = 1,at = c(0,(xlimChoice/2),xlimChoice),labels = c(0,(xlimChoice/2)/1000,xlimChoice/1000))
#   axis(side = 2,at = c(-.1, 0, .1, .2,.3), las = 2)
#   
#   for(j in 1:2){
#     
#     gType = c("E", "L")[j]
#     
#     TE.gr <- GRanges(seqnames = Rle(joinRep[[TEfam]]$genoName), 
#                      ranges = IRanges(start = joinRep[[TEfam]]$genoStart, end = joinRep[[TEfam]]$genoEnd))
#     
#     gapsOpen.gr  <- get(paste("gapsOpenCell", gType, ".gr", sep = ""))
#     
#     gapsOpen <- data.frame(chr = seqnames(gapsOpen.gr), start = start(gapsOpen.gr), end = end(gapsOpen.gr), Known = width(gapsOpen.gr))
#     gapsOpen[,TEfam] = countOverlaps(gapsOpen.gr, TE.gr)
#     
#     posStatBins = gapsOpen
#     
#     TEs_posStats <- covCalcPlot(lenChoice = lenChoice, repChoice = TEfam, repBins = posStatBins,repList = joinRep)
#     
#     rawCov <- (TEs_posStats$rawRepCov3 + TEs_posStats$rawRepCov5[(lenChoice+1):1])
#     bpFreq <- (TEs_posStats$baseFreq3 + TEs_posStats$baseFreq5[(lenChoice+1):1])
#     bpFreq[is.na(bpFreq)] <- 1
#     
#     rate <- rawCov/bpFreq
#     
#     # we can get both stats, 
#     cutSite <- unique(as.integer(10^(seq(0,as.integer(log10(lenChoice+1))+1,.05))))
#     cuted <- cut(1:(lenChoice+1), breaks = cutSite,right = F,ordered_result = T)
#     cutPos<- as.integer(apply(data.frame(cutSite[1:(length(cutSite)-1)], cutSite[2:(length(cutSite))] ), 1, mean))
#     Tab <- table(cuted)
#     
#     aggTEmean <- aggregate(rate, list(as.integer(cuted)), mean, simplify = F)
#     aggTEsd <- aggregate(rate, list(as.integer(cuted)), sd, simplify = F)
#     
#     aggBPfreq <- aggregate(bpFreq, list(as.integer(cuted)), sum, simplify = F)
#     aggBPraw <- aggregate(rawCov, list(as.integer(cuted)), sum, simplify = F)
#     
#     usePos <- cutPos[aggTEsd$Group.1]
#     
#     biasData <- get(paste(gType, TEfam, "Bias", sep =""))
#     SPbiasData <- smooth.spline((biasData$position[complete.cases(biasData)]),biasData$bias[complete.cases(biasData)],all.knots = TRUE)
#     biasPred <- predict(SPbiasData, usePos)
#     
#     sd <- (sqrt(bpFreq * TEs_posStats$p * (1 - TEs_posStats$p))/(bpFreq))
#     
#     shape.x <- c((usePos), rev((usePos)))
#     shape.y <- c((((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p) + (3 *sd[usePos]), 
#                  rev((((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p) - (3 * sd[usePos])))
#     
#     polygon(shape.x, shape.y, density =20,border = 0,col = j+1 ,angle = 135)
#     
#     lines((usePos), (((aggBPraw$x/aggBPfreq$x)/biasPred$y) - TEs_posStats$p), col = j + 1,lwd = 1 )
#     
#     abline(h = 0, lty = 2)   
#     
#   }
#   
# }
# 
# dev.off()
# 
# 









