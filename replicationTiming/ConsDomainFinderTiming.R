### so how can we look at shifts in gene density and new L1 density in domains 



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

# 
# web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
# con <- gzcon(url(web))
# txt <- readLines(con)
# refgene <- read.delim(textConnection(txt), header = FALSE)
# refgene.gr <- reduce(GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6])))
# 
# 
# rep <- rep_info(spec1 = spec1, genome = genome)
# 
# 
# joinRep <- list(old_L1 = rbind(rep$L1ME, rep$L1MD, rep$L1MC, rep$L1MB),
#                 new_L1 = rbind(rep$L1MA, rep$L1PB, rep$L1PA, rep$L1HS),
#                 Alu = rbind(rep$AluJ, rep$AluS, rep$AluY),
#                 Ancient = rbind(rep$MIR, rep$L2)
# )
# 
# #joinRep = rep
# ChromoChoice = "R"
# joinChromo <- c(list(`H1-hESC`[[ChromoChoice]]), list(HUVEC[[ChromoChoice]]), list(HepG2[[ChromoChoice]]), list(GM12878[[ChromoChoice]]), list(`HeLa-S3`[[ChromoChoice]]), list(K562[[ChromoChoice]]))
# names(joinChromo) <- c("H1hESC", "HUVEC", "HepG2", "GM12878", "HelaS3", "K562")
# names(joinChromo) <- paste(names(joinChromo), ChromoChoice, sep = "_")
# 
# chromoSample <- HUVEC
# 
# joinRepChromatin <- c(joinRep, chromoSample)



#### this is where we get our regions 
# 








waveList <- list.files("Data/UW_repliSeq/")
for(i in 1:length(waveList)){
  wave <- import(paste("Data/UW_repliSeq/", waveList[i], sep = ""), format = "bw")
  if(i == 1){
    waveSet <- data.frame(chr = as.character(seqnames(wave)), start = start(wave), end = end(wave))
  }
  waveName = gsub(pattern = "wgEncodeUwRepliSeq", replacement = "", waveList[i])
  waveName = gsub(pattern = "WaveSignalRep1.bigWig", replacement = "", waveName)
  waveSet[,waveName] = elementMetadata(wave)$score
}

waveSet[,4:ncol(waveSet)][waveSet[,4:ncol(waveSet)] < 6] <- NA
WaveSet <- scale(waveSet[,4:ncol(waveSet)])
waveMean <- rowMeans(WaveSet,na.rm = T)
waveSd <- apply(FUN = sd,X = WaveSet,MARGIN = 1,na.rm = T)

waveDomain <- rep(NA, nrow(WaveSet))
waveDomain[!is.na(waveMean)] <- 0
waveDomain[waveMean > 0 + (3 * waveSd)] <- 1
waveDomain[waveMean < 0 - (3 * waveSd)] <- -1
waveDomain <- data.frame(waveSet[,c("chr", "start", "end")], domain = waveDomain)

domainRanges <- NULL
chrName <- unique(as.character(waveDomain$chr))
for(chr in chrName){
  domains <- waveDomain[waveDomain$chr == chr,]
  domains <- domains[!is.na(domains$domain),]
  for(d in c(-1,0,1)){
    domainD <- domains[domains$domain == d,]
    DD <- domainD$start[2:nrow(domainD)] - domainD$end[1:(nrow(domainD)-1)]
    starts <- c(domainD[1,"start"] ,domainD[(1:nrow(domainD))[DD > 1] + 1,"start"])
    ends <- c(domainD[DD > 1,"end"], domainD[nrow(domainD),"end"])
    domainRanges = rbind(domainRanges, data.frame(chr = chr, start = starts, end = ends, domain = d))
  }
}

write.table(x = domainRanges, file = "Data/ConsTimingDomains", quote = F, sep = "\t", row.names = F, col.names = T)



slope <- NULL
for(i in 1001:100000){
  slope <- c(slope,(waveMean[i+100] - waveMean[i-100]) / 3)
}

plot(waveMean[10000:110000], type = "l")
par(new=T)
plot(x = 1001:100000, y = abs(slope), type = "l", xlim = c(0,100000))



high = 390000
low = 370000

plot(waveMean[low:high], type = "l", ylim = c(-2,2))
lines(waveMean[low:high] - (1 * waveSd[low:high]), col = 2)
lines(waveMean[low:high] + (1 * waveSd[low:high]), col = 2)
abline(h = .5)
abline(h=-.5)

par(new=T)
plot(waveDomain$domain[low:high],col = waveDomain$domain[low:high]+2, type = "b", ylim = c(-2,2))


plot(density(waveMean,na.rm = T))


### alternativly we could identify regions where the score doesn't change much 


# 
# 
# ERD <- domainRanges[domainRanges$domain == 1,]
# ERD.gr <- GRanges(seqnames = Rle(ERD$chr), 
#                   ranges = IRanges(start = ERD$start, end = ERD$end))
# 
# 
# 
# new_L1.gr <- GRanges(seqnames = Rle(joinRepChromatin$new_L1$genoName),
#                      ranges = IRanges(start = joinRepChromatin$new_L1$genoStart, end = joinRepChromatin$new_L1$genoEnd))
# 
# int.ref <- intersect(ERD.gr, refgene.gr)
# 
# int <- intersect(ERD.gr,new_L1.gr)
# 
# new_L1.OL <- as.matrix(findOverlaps(ERD.gr, int))
# new_L1.agg <- aggregate(x = width(int)[new_L1.OL[,2]], by = list(new_L1.OL[,1]), FUN = sum)
# 
# refgene.OL <- as.matrix(findOverlaps(ERD.gr, int.ref))
# refgene.agg <- aggregate(x = width(int.ref)[refgene.OL[,2]], by = list(refgene.OL[,1]), FUN = sum)
# 
# 
# 
# plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log( (new_L1.agg$x)/(width(ERD.gr)[new_L1.agg$Group.1])))
# plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log(new_L1.agg$x))
# abline(lm(log(new_L1.agg$x) ~ log(width(ERD.gr)[new_L1.agg$Group.1])))
# 
# ## so far the smaller the bin the higher the L1 density
# 
# plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log( (new_L1.agg$x)/(width(ERD.gr)[new_L1.agg$Group.1])))
# 
# plot(log(width(ERD.gr)[refgene.agg$Group.1]),log(refgene.agg$x))
# 
# 
# 
# A <- merge(new_L1.agg, refgene.agg, by = 1)
# plot(    ( (A$x.y)/(width(ERD.gr)[A$Group.1]))   ,   ( (A$x.x)/(width(ERD.gr)[A$Group.1])))
# cor(    ( (A$x.y)/(width(ERD.gr)[A$Group.1]))   ,   ( (A$x.x)/(width(ERD.gr)[A$Group.1])))
# 
# 
# plot( log10( (A$x.y))   ,  log10 ( (A$x.x)))
# 




head(waveDomain)


web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
chromLength <- read.delim(textConnection(txt), header = FALSE)




# lets set these up as a matrix

#rangesChr21 <- ranges[ranges[,1] == "chr21",]



bins <- binned.genome.reader(genome = "hg19", bin.size = c(1000), keep.rate = .9)
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
}


# so now we have a matrix and we can get a consensus on each state
binMat <- bins[,sampNames]


binMat <- as.matrix(binMat)
binMat[binMat == "ERD"] = "darkgreen"
binMat[binMat == "LRD"] = "blue"
binMat[binMat == "UTZ"] = "red"
binMat[binMat == "DTZ"] = "black"

states <- c("darkgreen", "blue", "red", "black")

Ymat <- matrix(data = 1:ncol(binMat), nrow = nrow(binMat), ncol = ncol(binMat),byrow = T)
#xlim = c(10000,13000)
high = 130000
low = 20000

plot(1,xlim=c(0,length(low:high)), ylim = c(0,16), col = binMat[,1], type = "n")

for(i in 1:ncol(binMat)){
  for(j in 1:length(states)){
    state1 <- rep(NA,(high - low)+1)
    state1[binMat[low:high,i] == states[j]] = i
    lines(state1, col = states[j], lwd = 20)
  }
 # points(Ymat[low:high,i], col = binMat[,i])
}


# maybe we should just take a consensus approach

binMat <- bins[,sampNames]

states <- c("ERD", "LRD", "UTZ", "DTZ")

stateDF <- matrix(data = 0, nrow = nrow(binMat), ncol = 5)
colnames(stateDF) <- c(states, "nonA")

for(i in 1:4){
  mat <- matrix(data = 0, nrow = nrow(binMat), ncol = ncol(binMat))
  mat[binMat == states[i]] <- 1
  stateDF[,states[i]] <- rowSums(mat)
}
mat <- matrix(data = 0, nrow = nrow(binMat), ncol = ncol(binMat))
mat[is.na(binMat)] <- 1
stateDF[,5] <- rowSums(mat)



A <- apply(stateDF, MARGIN = 1, FUN = max)


colFind <- function()



  colnames(stateDF)

colnames(stateDF[,stateDF[1,] == A[1]])



consensus = rep("", nrow(stateDF))
for(i in 1:5){
  B = stateDF[,i]
  consensus[B == A] = colnames(stateDF)[i]
}

conBin <- bins[1:4]
conBin$consensus = consensus
conBin$conScore <- A

# do we also want to add the mean score and sd across our timing set

# it would probably give us the most comprahensive information

## turn this into a range somehow by taking adjacent bins with same class


# show accuracy by counting how many equal bins we get between our consnesus and all others 

# and all others between each other 


# shoudl result in a heatmap 
# we could also include a probability paramater aswell, what degreen of beliefe do we place on 
# our assignment. 

conDomain <- conBin


domainRanges <- NULL
chrName <- unique(as.character(conDomain$chr))
for(chr in chrName){
  domains <- conDomain[conDomain$chr == chr,]
 # domains <- domains[!is.na(domains$domain),]
  for(d in c("nonA", "UTZ", "DTZ", "ERD", "LRD")){
    if(nrow(domains[domains$consensus == d,]) == 0){
      next
    }
    domainD <- domains[domains$consensus == d,]
    DD <- domainD$start[2:nrow(domainD)] - domainD$end[1:(nrow(domainD)-1)]
    starts <- c(domainD[1,"start"] ,domainD[(1:nrow(domainD))[DD > 1] + 1,"start"])
    ends <- c(domainD[DD > 1,"end"], domainD[nrow(domainD),"end"])
    domainRanges = rbind(domainRanges, data.frame(chr = chr, start = starts, end = ends, domain = d))
  }
}

chrName <- unique(as.character(conDomain$chr))
for(chr in chrName){
  domains <- domainRanges[domainRanges$chr == chr,]
  domains = domains[order(domains$start),]
  domainRanges[domainRanges$chr == chr,] = domains
}


write.table(x = domainRanges, file = "Data/ConsTimingDomains", quote = F, sep = "\t", row.names = F, col.names = T)


