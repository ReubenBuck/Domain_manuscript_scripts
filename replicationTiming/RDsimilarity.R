


# lets actually get the stuff from the constiuative domains script and look at percentage overlap of each domain class






rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")
#load("R_objects/chromStateCombined")


genome = "hg19"
spec1 = "Human"




bins <- binned.genome.reader(genome = "hg19", bin.size = c(1000), keep.rate = 0)
bins <- bins[[1]]

bins.gr <- GRanges(seqnames = Rle(bins$chr), 
                   ranges = IRanges(start = bins$start, end = bins$end))


timeFiles <- list.files(path = "Data/DNN_HMM_repliDomains/")
sampNames <- NULL
for( i in 1:length(timeFiles)){
  name <- paste(strsplit(timeFiles[i],"_")[[1]][3],strsplit(timeFiles[i],"_")[[1]][4] , sep = "_")
  ranges <- read.table(paste("Data/DNN_HMM_repliDomains/", timeFiles[i], sep =""))
  ranges.gr <- GRanges(seqnames = Rle(ranges[,1]),
                       ranges = IRanges(start = ranges[,2], end = ranges[,3]))
  ol <- as.matrix(findOverlaps(bins.gr, ranges.gr, type = c("within")))
  if(length(ol[duplicated(ol[,1]),]) != 0){
    print("duplicate overlap", i , name)
    break
  }
  bins[,name] <- NA
  class(bins[,name]) <- "character"
  bins[ol[,1],name] <- as.character(ranges[ol[,2],4])
  sampNames <- c(sampNames, name)
  
  colnames(ranges) <- c("chr", "start", "end", "domain")
  assign(x = name, value = ranges)
}



binMat <- bins[,sampNames]
binMat <- binMat[complete.cases(binMat),]
binMat[binMat == "UTZ"] <- "TTR"
binMat[binMat == "DTZ"] <- "TTR"

colSums(table(binMat[,1:2]))


# a high degree of overlap 

# do pairwise comparisons and get a percentage overlap for each region type 
# and do a fisher test

# and do corealtion of replication timing 

#sampleGrab <- expand.grid(data.frame(1:ncol(binMat), 1:ncol(binMat)))


TTRmat <- Fmat <- ERDmat <- LRDmat <- matrix(NA,nrow = ncol(binMat), 
                                             ncol = ncol(binMat), 
                                             dimnames = list(colnames(binMat),colnames(binMat))
                                             )
for(i in 1:ncol(binMat)){
  for(j in 1:ncol(binMat)){
    tab <- table(binMat[,c(i,j)])
    res <- diag(tab)/rowSums(tab)
    TTRmat[i,j] <- res["TTR"]
    ERDmat[i,j] <- res["ERD"]
    LRDmat[i,j] <- res["LRD"]
#    Fmat[i,j] <- fisher.test()
  }
}

colFun = colorRampPalette(c( "darkblue","lightblue", "white"))
image(t(matrix(1:20)),col = colFun(20), xaxt = "n", las = 2)
#image(ERDmat, col = colFun(20), zlim = c(0,1))


heatmap(TTRmat, scale = "none", zlim = c(0,1), col = colFun(20), margin = c(10,10))

heatmap(ERDmat, scale = "none", zlim = c(0,1), col = colFun(20), margin = c(10,10))

heatmap(LRDmat, scale = "none", zlim = c(0,1), col = colFun(20), margin = c(10,10))

# proportion overlap 

# what about 
# if we look at TTRs and how many samples are they shared across 

# if we find that most TTR regions are not unique 

numMat <- sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
numMat[binMat == "TTR"] = 1
CS <- colSums(numMat)

RS <- rowSums(numMat)

sampRS <- NULL
sampMats <- NULL
for(z in 1:10){
sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
for(j in 1:ncol(sampMat)){
  A1 <- (1:nrow(numMat))[numMat[,j] == 1]
  A2 <- (1:length(A1))[A1[2:length(A1)] - A1[1:(length(A1)-1)]>1]
  A3 <- A2 - c(0,A2[1:(length(A2)-1)])
  
  S1 <- sample(x = A3, size = as.integer(nrow(numMat)/mean(A3)), replace = TRUE)
  S2 = NULL
  for(i in 1:length(S1)){S2 <- c(S2,sum(S1[1:i]))}
  S2 <- S2[S2 < nrow(sampMat)]
  S3 <- sample(x = 1:length(S2),length(A3),replace = F)
  
  intervals <- NULL
  for(i in 1:length(S3)){
    intervals <- c(intervals,(S2[S3] + 1 - S1[S3])[i]:(S2[S3] + 1)[i])
  }
  sampMat[intervals,j] = 1
}

colnames(sampMat) <- colnames(binMat)
sampMats <- c(sampMats, list(sampMat))
sampRS <- rbind(sampRS, rowSums(sampMat))
}



pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/replicationDomain/TTRpercentOL.pdf", 
    height = 5, width =5)
h1 <- hist(RS[RS > 0], breaks = seq(.5,16.5), plot = FALSE)
h2 <- hist(sampRS[1,][sampRS[1,] > 0], breaks = seq(.5,16.5), plot = FALSE)

plot(h2$mids,h2$counts/sum(sampRS[1,][sampRS[1,] > 0]), type = "n", lwd = 3, pch = 16, xlab = "number of samples", ylab = "genome domain fraction")
for(i in 1:10){
  h2 <- hist(sampRS[i,][sampRS[i,] > 0], breaks = seq(.5,16.5), plot = FALSE)
  lines(h2$mids,h2$counts/sum(sampRS[i,][sampRS[i,] > 0]), type = "l", col = 8, lwd = 1)
}
points(h1$mids,h1$counts/sum(RS[RS>0]), type = "b", lwd = 3, pch = 16)
legend("topright", legend = c("observed", "expected"), fill = c(1,8), bty = "n")
dev.off()
# this tells us the level of expected overlap


# we could look at the level of expected between two pairwise samples and get a measure on variance 

### lets get some dinner and finalise these plots where we establish decent overlap between regions. 
### is it worth doing pairwise overlap across our samples?



numMat <- sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
numMat[binMat == "ERD"] = 1
CS <- colSums(numMat)

RS <- rowSums(numMat)

sampRS <- NULL
sampMats <- NULL
for(z in 1:10){
  sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
  for(j in 1:ncol(sampMat)){
    A1 <- (1:nrow(numMat))[numMat[,j] == 1]
    A2 <- (1:length(A1))[A1[2:length(A1)] - A1[1:(length(A1)-1)]>1]
    A3 <- A2 - c(0,A2[1:(length(A2)-1)])
    
    S1 <- sample(x = A3, size = as.integer(nrow(numMat)/mean(A3)), replace = TRUE)
    S2 = NULL
    for(i in 1:length(S1)){S2 <- c(S2,sum(S1[1:i]))}
    S2 <- S2[S2 < nrow(sampMat)]
    S3 <- sample(x = 1:length(S2),length(A3),replace = F)
    
    intervals <- NULL
    for(i in 1:length(S3)){
      intervals <- c(intervals,(S2[S3] + 1 - S1[S3])[i]:(S2[S3] + 1)[i])
    }
    sampMat[intervals,j] = 1
  }
  
  colnames(sampMat) <- colnames(binMat)
  sampMats <- c(sampMats, list(sampMat))
  sampRS <- rbind(sampRS, rowSums(sampMat))
}



pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/replicationDomain/ERDpercentOL.pdf", 
    height = 5, width =5)
h1 <- hist(RS[RS > 0], breaks = seq(.5,16.5), plot = FALSE)
h2 <- hist(sampRS[1,][sampRS[1,] > 0], breaks = seq(.5,16.5), plot = FALSE)

plot(h2$mids,h2$counts/sum(sampRS[1,][sampRS[1,] > 0]), type = "n", lwd = 3, pch = 16, xlab = "number of samples", ylab = "genome domain fraction")
for(i in 1:10){
  h2 <- hist(sampRS[i,][sampRS[i,] > 0], breaks = seq(.5,16.5), plot = FALSE)
  lines(h2$mids,h2$counts/sum(sampRS[i,][sampRS[i,] > 0]), type = "l", col = 8, lwd = 1)
}
points(h1$mids,h1$counts/sum(RS[RS>0]), type = "b", lwd = 3, pch = 16)
legend("topright", legend = c("observed", "expected"), fill = c(1,8), bty = "n")
dev.off()
# this tells us the level of expected overlap



#### LRD 




numMat <- sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
numMat[binMat == "LRD"] = 1
CS <- colSums(numMat)

RS <- rowSums(numMat)

sampRS <- NULL
sampMats <- NULL
for(z in 1:10){
  sampMat <- matrix(0,nrow = nrow(binMat), ncol = ncol(binMat))
  for(j in 1:ncol(sampMat)){
    A1 <- (1:nrow(numMat))[numMat[,j] == 1]
    A2 <- (1:length(A1))[A1[2:length(A1)] - A1[1:(length(A1)-1)]>1]
    A3 <- A2 - c(0,A2[1:(length(A2)-1)])
    
    S1 <- sample(x = A3, size = as.integer(nrow(numMat)/mean(A3)), replace = TRUE)
    S2 = NULL
    for(i in 1:length(S1)){S2 <- c(S2,sum(S1[1:i]))}
    S2 <- S2[S2 < nrow(sampMat)]
    S3 <- sample(x = 1:length(S2),length(A3),replace = F)
    
    intervals <- NULL
    for(i in 1:length(S3)){
      intervals <- c(intervals,(S2[S3] + 1 - S1[S3])[i]:(S2[S3] + 1)[i])
    }
    sampMat[intervals,j] = 1
  }
  
  colnames(sampMat) <- colnames(binMat)
  sampMats <- c(sampMats, list(sampMat))
  sampRS <- rbind(sampRS, rowSums(sampMat))
}



pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/replicationDomain/LRDpercentOL.pdf", 
    height = 5, width =5)
h1 <- hist(RS[RS > 0], breaks = seq(.5,16.5), plot = FALSE)
h2 <- hist(sampRS[1,][sampRS[1,] > 0], breaks = seq(.5,16.5), plot = FALSE)

plot(h2$mids,h2$counts/sum(sampRS[1,][sampRS[1,] > 0]), type = "n", lwd = 3, pch = 16, xlab = "number of samples", ylab = "genome domain fraction")
for(i in 1:10){
  h2 <- hist(sampRS[i,][sampRS[i,] > 0], breaks = seq(.5,16.5), plot = FALSE)
  lines(h2$mids,h2$counts/sum(sampRS[i,][sampRS[i,] > 0]), type = "l", col = 8, lwd = 1)
}
points(h1$mids,h1$counts/sum(RS[RS>0]), type = "b", lwd = 3, pch = 16)
legend("topright", legend = c("observed", "expected"), fill = c(1,8), bty = "n")
dev.off()
# this tells us the level of expected overlap

