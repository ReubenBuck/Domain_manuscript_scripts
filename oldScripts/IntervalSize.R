###### If i cant get the function to work propoerly i can at least re write some of the stuff here


## it will be cool to compare to intergenic regions 




rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
source("Domain_manuscript_scripts/functions.R")

library(GenomicRanges)
library(rtracklayer)


spec1 <- "Human"
genome = "hg19"


web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)
colnames(refgene) <- c("bin", "ID", "chr", "strand", "start", "end", "txnStart", "txnEnd", "exonCount",
                       "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")


tsd <- read.table(file = "~/Desktop/quilt/tsdFinder/hg19_chrKnown.tsd.txt", skip = 1)
colnames(tsd) <- c("chr", "start", "end", "leftStart", 
                   "rightStart", "tsdLength", "repFam", 
                   "repType", "repStart", "repEnd", "repLeft",
                   "repFam2", "repType2" , "tsdSeq", "strand")


intronKeep.gr <- filterIntron(refgene)
refgene_gap.gr <- filterIntergenic(refgene)

# if we can start comparing the differnt TE classes that would be cool


old_L1name <- c("L1ME", "L1MD", "L1MC", "L1MB")
new_L1name <- c("L1MA", "L1PB", "L1PA", "L1HS")
Ancientname <- c("L2", "MIR")
Aluname <- c("AluJ", "AluS", "AluY")
TEname <- c("old_L1", "new_L1", "Ancient", "Alu")

for(i in 1:length(TEname)){
  FAMnames <- get(paste(TEname[i], "name", sep = ""))
  print(FAMnames)
  TEgroup = NULL
  for(j in 1:length(FAMnames)){
    TEgroup <- rbind(TEgroup,tsd[grep(FAMnames[j], tsd$repFam),])
    print(FAMnames[j])
  }
  assign(x = TEname[i], value = TEgroup)
}



IntronSamp <- sample(width(intronKeep.gr), 100000,prob = width(intronKeep.gr), replace = T)
IntronHist <- hist(log2(IntronSamp), breaks=seq(0,25,.1))

GapSamp <- sample(width(refgene_gap.gr), 100000,prob = width(refgene_gap.gr), replace = T)
GapHist <- hist(log2(GapSamp), breaks = seq(0,25,.1))


#plot(c(0,30), c(0,.06), type = "n")


#plot(c(0,(1000000)), c(0,.06), type = "n")
pdf(file = "plots/TSD_intervalSize.pdf", width = 10,height = 20)
layout(matrix(c(1,2,3,4,5,6,7,8),nrow = 4))


minTSD <- c(6,9,6,9)
min3prime <-c(1,6000,1,280)
cols <- c("red", "purple", "darkblue","darkgreen")
for(i in 1:length(TEname)){
  TEFAM <- get(TEname[i])
  
  TEgroupTSD <- TEFAM[TEFAM[,"tsdLength"] >= minTSD[i] & TEFAM[,"repEnd"] > min3prime[i],]
  tsd.gr <- GRanges(seqnames = Rle(TEgroupTSD$chr),
                      ranges = IRanges(start = TEgroupTSD$start, end = TEgroupTSD$end)
  )
  
  fOLtsdINTRON <- as.matrix(findOverlaps(tsd.gr, intronKeep.gr))
  myhistINTRON <- hist(log2(width(intronKeep.gr)[fOLtsdINTRON[,2]]), breaks = IntronHist$breaks,plot = F)
  
  x = (myhistINTRON$mids)
  y = (myhistINTRON$counts/IntronHist$counts)
  df <- data.frame(x = x, y = y)
  df$y[is.infinite(df$y)] <- NA
  df$y =df$y   #/sum(df$y,na.rm = T)
  
  n = sum(myhistINTRON$counts, na.rm = T)
  p = (IntronHist$counts[IntronHist$counts>0])/sum((IntronHist$counts[IntronHist$counts>0]))
  
  
  df$tsdMean[IntronHist$counts>0] = (n*p)/IntronHist$counts[IntronHist$counts>0]
  df$sd[IntronHist$counts>0] <- sqrt(n*p * (1-p))/IntronHist$counts[IntronHist$counts>0]
  
  
  df = df[complete.cases(df),]
  lo <- loess(df$y ~ df$x)
  plot(df$x, log2((df$y) * (2^df$x) / sum((df$y) * (2^df$x))), col = cols[i], pch = 16, lwd = 3,xaxt= "n", yaxt = "n",
       xlab = "interval size (kb)", ylab = "adjusted insertion events (%)", main = paste(TEname[i], "Intron"), cex = .5, 
       xlim = c(8,25), ylim = c(-20,-2))
  
  axis(side=1, at=seq(0,100,2), labels = round((2^seq(0,100,2))/1000, digits = 2))
  axis(side=2, at=seq(-100,100,5), labels = round((2^seq(-100,100,5))*100, digits = 4) )
#  lines(df$x, (df$tsdMean * (2^df$x) / sum((df$y) * (2^df$x))), col = 1, lwd = 3, lty = 1)
  lines(df$x, log2((df$tsdMean + (2*df$sd)) * (2^df$x) / sum((df$y) * (2^df$x)))  , col = 1, lwd = 1, lty = 1)
  lines(df$x, log2((df$tsdMean - (2*df$sd)) * (2^df$x) / sum((df$y) * (2^df$x)))  , col = 1, lwd = 1, lty = 1)
  
  
}


for(i in 1:length(TEname)){
  TEFAM <- get(TEname[i])
  
  TEgroupTSD <- TEFAM[TEFAM[,"tsdLength"] >= minTSD[i] & TEFAM[,"repEnd"] > min3prime[i],]
  tsd.gr <- GRanges(seqnames = Rle(TEgroupTSD$chr),
                    ranges = IRanges(start = TEgroupTSD$start, end = TEgroupTSD$end)
  )
  
  fOLtsdINTERGENIC <- as.matrix(findOverlaps(tsd.gr, refgene_gap.gr))
  myhistINTERGENIC <- hist(log2(width(refgene_gap.gr)[fOLtsdINTERGENIC[,2]]), breaks = GapHist$breaks,plot = F)
  
  x = (myhistINTERGENIC$mids)
  y = (myhistINTERGENIC$counts/GapHist$counts)
  df <- data.frame(x = x, y = y)
  df$y[is.infinite(df$y)] <- NA
  df$y =df$y   #/sum(df$y,na.rm = T)
  
  n = sum(myhistINTERGENIC$counts, na.rm = T)
  p = (GapHist$counts[GapHist$counts>0])/sum((GapHist$counts[GapHist$counts>0]))
  
  
  df$tsdMean[GapHist$counts>0] = (n*p)/GapHist$counts[GapHist$counts>0]
  df$sd[GapHist$counts>0] <- sqrt(n*p * (1-p))/GapHist$counts[GapHist$counts>0]
  
  
  df = df[complete.cases(df),]
  lo <- loess(df$y ~ df$x)
  plot(df$x, log2(((df$y) * (2^df$x)) / sum((df$y) * (2^df$x) )), col = cols[i], pch = 16, lwd = 3,xaxt= "n", yaxt = "n",
       xlab = "interval size (kb)", ylab = "adjusted insertion events (%)", main = paste(TEname[i], "intergenic"), cex = .5,
       xlim = c(8,25), ylim = c(-20,-2))
  
  axis(side=1, at=seq(0,100,2), labels = round((2^seq(0,100,2))/1000, digits = 2))
  axis(side=2, at=seq(-100,100,5), labels = round((2^seq(-100,100,5))*100, digits = 4) )
  #lines(df$x, (df$tsdMean * (2^df$x) / sum((df$y) * (2^df$x))), col = 1, lwd = 3, lty = 1)
  lines(df$x, log2((df$tsdMean + (2*df$sd)) * (2^df$x) / sum((df$y) * (2^df$x)))  , col = 1, lwd = 1, lty = 1)
  lines(df$x, log2((df$tsdMean - (2*df$sd)) * (2^df$x) / sum((df$y) * (2^df$x)))  , col = 1, lwd = 1, lty = 1)
  
  #lines(df$x, log2((df$tsdMean - (2*df$sd))  * (2^df$x))  , col = 1, lwd = 3, lty = 3)
}


# another thing will be to see how intron lengths associate with 


dev.off()






pdf(file = "plots/TSD_intervalSizeFlat.pdf", width = 10,height = 20)
layout(matrix(c(5,6,7,8,1,2,3,4),nrow = 4))


minTSD <- c(6,9,6,9)
min3prime <-c(1,6000,1,280)
cols <- c("red", "purple", "darkblue","darkgreen")
for(i in 1:length(TEname)){
  TEFAM <- get(TEname[i])
  
  TEgroupTSD <- TEFAM[TEFAM[,"tsdLength"] >= minTSD[i] & TEFAM[,"repEnd"] > min3prime[i],]
  tsd.gr <- GRanges(seqnames = Rle(TEgroupTSD$chr),
                    ranges = IRanges(start = TEgroupTSD$start, end = TEgroupTSD$end)
  )
  
  fOLtsdINTRON <- as.matrix(findOverlaps(tsd.gr, intronKeep.gr))
  myhistINTRON <- hist(log2(width(intronKeep.gr)[fOLtsdINTRON[,2]]), breaks = IntronHist$breaks,plot = F)
  
  x = (myhistINTRON$mids)
  y = (myhistINTRON$counts/IntronHist$counts)
  df <- data.frame(x = x, y = y)
  df$y[is.infinite(df$y)] <- NA
  df$y =df$y   #/sum(df$y,na.rm = T)
  n = sum(myhistINTRON$counts, na.rm = T)
  p = (IntronHist$counts[IntronHist$counts>0])/sum((IntronHist$counts[IntronHist$counts>0]))
  
  
  df$tsdMean[IntronHist$counts>0] = (n*p)/IntronHist$counts[IntronHist$counts>0]
  df$sd[IntronHist$counts>0] <- sqrt(n*p * (1-p))/IntronHist$counts[IntronHist$counts>0]
  
  
  df = df[complete.cases(df),]
  lo <- loess(df$y ~ df$x)
  plot(df$x, (df$y) / sum(df$y), col = cols[i], pch = 16, lwd = 3,xaxt= "n",
       xlab = "interval size (kb)", ylab = "insertion events per available space (%)", main = paste(TEname[i], "Intron"), cex = .5,
       ylim = c(0,.025), xlim = c(8,25)
       )
  
  axis(side=1, at=seq(0,100,2), labels = round((2^seq(0,100,2))/1000, digits = 2))
#  axis(side=2, at=seq(-100,100,5), labels = round((2^seq(-100,100,5))*100, digits = 4) )
  lines(df$x, (df$tsdMean)/sum((df$y)), col = 1, lwd = 3, lty = 1)
  lines(df$x, (df$tsdMean + (2*df$sd))  / sum((df$y))  , col = 1, lwd = 1, lty = 1)
  lines(df$x, (df$tsdMean - (2*df$sd))  / sum((df$y))  , col = 1, lwd = 1, lty = 1)
  
  pred <- predict(lo,df$x, se = T)
  lines(df$x,  pred$fit / sum(df$y), col = cols[i], lwd = 3, lty = 1)
  polygon(c(df$x, rev(df$x)), c((pred$fit + (2 * pred$se.fit) )/ sum(df$y), rev((pred$fit - (2 * pred$se.fit) )/ sum(df$y))),
          col = cols[i],density = 20, border = NA)
  
  
}


for(i in 1:length(TEname)){
  
  
  TEFAM <- get(TEname[i])
  
  TEgroupTSD <- TEFAM[TEFAM[,"tsdLength"] >= minTSD[i] & TEFAM[,"repEnd"] > min3prime[i],]
  tsd.gr <- GRanges(seqnames = Rle(TEgroupTSD$chr),
                    ranges = IRanges(start = TEgroupTSD$start, end = TEgroupTSD$end)
  )
  
  fOLtsdINTERGENIC <- as.matrix(findOverlaps(tsd.gr, refgene_gap.gr))
  myhistINTERGENIC <- hist(log2(width(refgene_gap.gr)[fOLtsdINTERGENIC[,2]]), breaks = GapHist$breaks,plot = F)
  
  
  x = (myhistINTERGENIC$mids)
  y = (myhistINTERGENIC$counts/GapHist$counts)
  df <- data.frame(x = x, y = y)
  df$y[is.infinite(df$y)] <- NA
  df$y = df$y   #/sum(df$y,na.rm = T)
  
  n = sum(myhistINTERGENIC$counts, na.rm = T)
  p = (GapHist$counts[GapHist$counts>0])/sum((GapHist$counts[GapHist$counts>0]))
  
  
  df$tsdMean[GapHist$counts>0] = (n*p)/GapHist$counts[GapHist$counts>0]
  df$sd[GapHist$counts>0] <- sqrt(n*p * (1-p))/GapHist$counts[GapHist$counts>0]
  
  
  df = df[complete.cases(df),]
  lo <- loess(df$y ~ df$x)
  plot(df$x, (df$y) / sum(df$y), col = cols[i], pch = 16, lwd = 3,xaxt= "n",
       xlab = "interval size (kb)", ylab = "adjusted insertion events (%)", main = paste(TEname[i], "intergenic"), cex = .5,
       xlim = c(8,25))
  
  axis(side=1, at=seq(0,100,2), labels = round((2^seq(0,100,2))/1000, digits = 2))
 # axis(side=2, at=seq(-100,100,5), labels = round((2^seq(-100,100,5))*100, digits = 4) )
  lines(df$x, df$tsdMean  / sum((df$y)), col = 1, lwd = 3, lty = 1)
  lines(df$x, (df$tsdMean + (2*df$sd))  / sum((df$y))   , col = 1, lwd = 1, lty = 1)
  lines(df$x, (df$tsdMean - (2*df$sd))  / sum((df$y))   , col = 1, lwd = 1, lty = 1)
  
  pred <- predict(lo,df$x, se = T)
  lines(df$x,  pred$fit / sum(df$y), col = cols[i], lwd = 3, lty = 1)
  polygon(c(df$x, rev(df$x)), c((pred$fit + (2 * pred$se.fit) )/ sum(df$y), rev((pred$fit - (2 * pred$se.fit) )/ sum(df$y))),
          col = cols[i],density = 20, border = NA)
  
}


dev.off()




plot(GapHist)
# we have to ask how many bases are 
2^mean(log2(width(refgene_gap.gr)))
#hist(GapSamp)
# for every base we get all the intervals that are bigger
# so do we do it per base pair or per interval of size x

# at this point im thinking per interval at size x
sizesMeanINTERGENIC <- NULL
for(i in seq(0,80000,100)){
  sizesMeanINTERGENIC <- c(sizesMeanINTERGENIC, 2^mean(log2(width(refgene_gap.gr)[width(refgene_gap.gr) > i])))
}


plot((seq(0,80000,100)/2),(predict(lo,newdata = log2(sizesMeanINTERGENIC))/ sum(df$y)) /  (df$tsdMean[1]  / sum((df$y))),
     main = paste(TEname[i], "intergenic"), xlab = "position (bp)", ylab = "interval size bias")

# we need to work out how to scale this 

# pased on size bias at position one 


sizesMeanINTRON <- NULL
for(i in seq(0,40000,100)){
  sizesMeanINTRON <- c(sizesMeanINTRON, 2^mean(log2(width(intronKeep.gr)[width(intronKeep.gr) > i])))
}
 

plot((seq(0,40000,100)/2),(predict(lo,newdata = log2(sizesMeanINTRON))/ sum(df$y)) /  (df$tsdMean[1]  / sum((df$y))),
     main = paste(TEname[i], "intron"), xlab = "position (bp)", ylab = "interval size bias")




plot(GapHist)
log2(10000)


# I guess half the inteergenic gneome is within 10000 bp of a gene
# an interesting stat 

# here we calculate the biasing factors at each position 
# maybe rather than sampling TSDs we sample repeat bases 
# how many bp are in ranges of a certain size 


