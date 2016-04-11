
#### here we can look at replication boundrreies identified by the human analysis in bioinformatics
### we can compare it with the big wig files

# this will give us a sense in human of which level 

rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")
load("R_objects/chromStateCombined")


genome = "hg19"
spec1 = "Human"



rep <- rep_info(spec1 = spec1, genome = genome)


joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
                new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
                Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
                Ancient = rbind(rep$MIR, rep$L2)
)

 #joinRep = rep
ChromoChoice = "R"
joinChromo <- c(list(`H1-hESC`[[ChromoChoice]]), list(HUVEC[[ChromoChoice]]), list(HepG2[[ChromoChoice]]), list(GM12878[[ChromoChoice]]), list(`HeLa-S3`[[ChromoChoice]]), list(K562[[ChromoChoice]]))
names(joinChromo) <- c("H1hESC", "HUVEC", "HepG2", "GM12878", "HelaS3", "K562")
names(joinChromo) <- paste(names(joinChromo), ChromoChoice, sep = "_")

chromoSample <- HUVEC

joinRepChromatin <- c(joinRep, chromoSample)





#### this is where we get our regions 
# 
# 
# fileNames <- list.files("Data/DNN_HMM_repliDomains/")
# cellName <- as.matrix(as.data.frame(strsplit(fileNames, "_"))[3,])[1:length(fileNames)]
# RTdomain <- NULL
# for(i in 1:length(fileNames)){
#   cell <- read.table(paste("Data/DNN_HMM_repliDomains/", fileNames[i], sep = ""),header = F,
#                      colClasses = c("character", "integer", "integer", "character"))
#   colnames(cell) <- c("chr", "start", "end", "domain")
#   RTdomain <- c(RTdomain,list(cell))
# }
# names(RTdomain) <- cellName
# 
# # really we want to pull out a boundry line?
# # basicly the start of a down domain and the end of an up domaintwisted around
# cellType <- "Huvec"
# 
TTR <- RTdomain[[cellType]][RTdomain[[cellType]]$domain %in% c("UTZ","DTZ") ,]
Edomain <- RTdomain[[cellType]][RTdomain[[cellType]]$domain == "ERD" ,]
Ldomain <- RTdomain[[cellType]][RTdomain[[cellType]]$domain == "LRD" ,]


RTdomains <- read.table("Data/ConsTimingDomains", header = T)
RTdomains <- RTdomains[complete.cases(RTdomains),]
Edomain <- RTdomains[RTdomains$domain == "ERD",]
Ldomain <- RTdomains[RTdomains$domain == "LRD",]



boxplot(list(TTR = TTR$end - TTR$start, Edomain = Edomain$end - Edomain$start,Ldomain = Ldomain$end- Ldomain$start),log = "y")


TTRchr <- c(TTR$chr[TTR$domain == "UTZ"], TTR$chr[TTR$domain == "DTZ"])
TTReb <- c(TTR$end[TTR$domain == "UTZ"], TTR$start[TTR$domain == "DTZ"])
TTRlb <- c(TTR$start[TTR$domain == "UTZ"], TTR$end[TTR$domain == "DTZ"])

TTRrotate <- rep(2, nrow(TTR))
TTRrotate[1:nrow(TTR[TTR$domain == "UTZ",])] = 1


binSize <- 50000
regionSize = 1000000

# GenomeBins <- binned.genome.reader(genome = genome,bin.size = binSize,keep.rate = .9)[[1]]
# BinnedReps <- binSort(repList = joinRepChromatin, bins = GenomeBins,TE.names = names(joinRepChromatin),repType = c(rep("repeats", length(joinRep)), rep("chromatin", length(chromoSample))))
# sampledMeans <- NULL
# for(i in 1:10000){
#   GenomeSample <- BinnedReps$rates[sample(x = 1:nrow(BinnedReps$rates), size = nrow(TTR),replace = F),]
#   sampledMeans <- rbind(sampledMeans,colMeans(GenomeSample[,5:(length(joinRepChromatin) + 4)]))
# }
# 
# 
# repMean <- colMeans(sampledMeans)
# repSD <- apply(sampledMeans, 2,sd)


# the ecdf can give us estimates of the number in whihch 5% of our data are above 


# so we can fit paramaters here to actually get a curve representative of the porbability distribution
# or we could do a data driven aproach and look at the actuall probabilitis


earlyAnno <- boundryAnotation(repList = joinRepChromatin, binSize = binSize, regionSize = regionSize, boundryLine = TTReb, 
                 repTypes = c(rep("repeats", length(joinRep)), rep("chromatin", length(chromoSample))),
                 boundryChr = TTRchr)
 

earlyAnno <- boundryAnotation(repList = joinRepChromatin[1:4], binSize = binSize, regionSize = regionSize, boundryLine = TTReb, 
                              repTypes = c(rep("repeats", length(joinRep))),
                              boundryChr = TTRchr)

# sampledMeans <- NULL
# for(i in 1:10000){
#   GenomeSample <- sample(x = earlyAnno$Ancient, size = nrow(TTR),replace = F)
#   sampledMeans <- c(sampledMeans,mean(GenomeSample))
# }


               
lateAnno <- boundryAnotation(repList = joinRepChromatin, binSize = binSize, regionSize = regionSize, boundryLine = TTRlb, 
                              repTypes = c(rep("repeats", length(joinRep)), rep("chromatin", length(chromoSample))),
                              boundryChr = TTRchr)

lateAnno <- boundryAnotation(repList = joinRepChromatin[1:4], binSize = binSize, regionSize = regionSize, boundryLine = TTRlb, 
                             repTypes = c(rep("repeats", length(joinRep)), rep("chromatin", length(chromoSample)))[1:4],
                             boundryChr = TTRchr)

#### wave file 


whichEB <- GRanges(seqnames=Rle(TTRchr), 
                   ranges=IRanges(start = TTReb - (regionSize/2), end = TTReb + (regionSize/2))
)
bigEB <- import("Data/UW_repliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format = "bw", which = whichEB)

whichLB <- GRanges(seqnames=Rle(TTRchr), 
                   ranges=IRanges(start = TTRlb - (regionSize/2), end = TTRlb + (regionSize/2))
)
bigLB <- import("Data/UW_repliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format = "bw", which = whichLB)



waveEarly <- boundryScoreAnotation(repList= list(bigEB),TE.names = "wave", repTypes = "GRange",metadata = "score",binSize = binSize,regionSize = regionSize,boundryLine = TTReb,boundryChr = TTRchr)
waveLate <- boundryScoreAnotation(repList = list(bigLB),TE.names = "wave", repTypes = "GRange",metadata = "score",binSize = binSize,regionSize = regionSize,boundryLine = TTRlb,boundryChr = TTRchr)
waveEarlyRot = waveEarly[[1]]
waveLateRot = waveLate[[1]]
waveEarlyRot[TTRrotate == 1,] <- rev(waveEarlyRot[TTRrotate == 1,])
waveLateRot[TTRrotate == 1,] <- rev(waveLateRot[TTRrotate == 1,])
WaveList <- list(waveEarlyRot, waveLateRot)

#### plotting things from beofre

pdf(file = "plots/repliBoundry/family_level_HUVEC_TEs_and_Segments_new_stat2.pdf", onefile = TRUE, width = 15, height = (15/4))
for(TEs in names(joinRepChromatin)){
  TE_choice <- TEs
  matList <- c(list(earlyAnno[[TE_choice]]/10000),list(lateAnno[[TE_choice]]/10000))
  
  
  layout(matrix(1:(2), nrow = 1))
  rotateTEmat <- NULL
  for(i in 1:2){
    TEmat = matList[[i]]
    #TEmat <- t(apply(X = matList[[i]], MARGIN = 1, smooth))
    TEmatUnflip <- TEmat
    TEmat[TTRrotate == 1,] <- rev(TEmat[TTRrotate ==1,])
    
    clust_pattern <- hclust(dist(TEmat))
    TEmat <- TEmat[clust_pattern$order,]
    rotateTEmat <- c(rotateTEmat, list(TEmat))
  }
  
  for(i in 1:2){
    
    sampledMeans <- NULL
    for(s in 1:10000){
      GenomeSample <- sample(x = matList[[i]], size = nrow(TTR),replace = F)
      sampledMeans <- c(sampledMeans,mean(GenomeSample))
    }
    repMean <- mean(sampledMeans)
    repSD <- sd(sampledMeans)
    
    
    y = colMeans(rotateTEmat[[i]])
    yMax <- max(c(max(abs(repMean - colMeans(rotateTEmat[[1]])) + repMean),
                  max(abs(repMean - colMeans(rotateTEmat[[1]])) + repMean), 
                  repMean + 5*repSD))
    waveMax <- max(c(colMeans(WaveList[[1]]), colMeans(WaveList[[2]])))
    waveMin <- min(c(colMeans(WaveList[[1]]), colMeans(WaveList[[2]])))
    
    par(mar=c(5,5,5,5))
    plot(y , type = "l", xaxt = "n", xlab= "", ylab = paste(TE_choice),col = 1 , ylim = c(repMean - (yMax - repMean),yMax), lwd = 3)
    abline(v=ncol(TEmat)/2 + .5)
    abline(h = repMean)
    abline(h = repMean + (2*repSD), lty = 2)
    abline(h = repMean - (2*repSD), lty = 2)
    
    par(new=T)
    plot(colMeans(WaveList[[i]]), xaxt  = "n", yaxt = "n", ylim = c(waveMin,waveMax), type = "l", lwd = 3, col = 2, ylab = "", xlab = "")
    axis(4, at = c(waveMin, waveMax), labels = c("Late", "Early"), col = 2)
    mtext("replication timing", side = 4, line = 2.5,col = 2)
    
  #  par(mar=c(5,5,0,5))
  #  image(t(rotateTEmat[[i]]), col = grey(seq(1,0,-.005)), xaxt = "n", yaxt = "n")
  #  abline(v = .5)
  #  axis(side = 1,at = c(0,.5,1), labels = c(paste(-regionSize/2/1000, "kb"), c("early boundry", "late boundry")[i], paste(regionSize/2/1000, "kb")))
  }
  
}

plot(y , type = "n",xaxt = "n" ,xlab= "", ylab = paste(TE_choice, c("Early", "Late")[i]),col = 1 , ylim = c(repMean - (yMax - repMean),yMax), lwd = 3)
axis(side = 1, at = c(0,10,20,30,40), c(-1000, -500, 0, 500, 1000))
plot(y , type = "n",xaxt = "n" ,xlab= "", ylab = paste(TE_choice, c("Early", "Late")[i]),col = 1 , ylim = c(repMean - (yMax - repMean),yMax), lwd = 3)
axis(side = 1, at = c(0,10,20,30,40), c(-1000, -500, 0, 500, 1000))


dev.off()
### will be interesting to do the late boarder aswell! 


### both borders side by side


#### wave signal

# we extract a section 

whichEB <- GRanges(seqnames=Rle(TTRchr), 
                 ranges=IRanges(start = TTReb - (regionSize/2), end = TTReb + (regionSize/2))
)
bigEB <- import("Data/UW_repliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format = "bw", which = whichEB)

whichLB <- GRanges(seqnames=Rle(TTRchr), 
                   ranges=IRanges(start = TTRlb - (regionSize/2), end = TTRlb + (regionSize/2))
)
bigLB <- import("Data/UW_repliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format = "bw", which = whichLB)




#### GC curve



A <- boundryScoreAnotation(repList = list(bigEB),TE.names = "wave", repTypes = "GRange",metadata = "score",binSize = 50000,regionSize = 2000000,boundryLine = TTReb,boundryChr = TTRchr)
Alate <- boundryScoreAnotation(repList = list(bigLB),TE.names = "wave", repTypes = "GRange",metadata = "score",binSize = 50000,regionSize = 2000000,boundryLine = TTRlb,boundryChr = TTRchr)

layout(c(1,2))
plot(c(colMeans(A[[1]]),colMeans(Alate[[1]])), type = "l")
plot(colMeans(Alate[[1]]), col = 2)






