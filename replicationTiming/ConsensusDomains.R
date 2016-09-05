
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


colnames(stateDF[,stateDF[1,] == A[1]])



consensus = rep("", nrow(stateDF))
for(i in 1:5){
  B = stateDF[,i]
  consensus[B == A] = colnames(stateDF)[i]
}

conBin <- bins[1:4]
conBin$consensus = consensus
conBin$conScore <- A


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
    if(nrow(domainD[domainD$conScore >= 12,]) == 0){
      next
    }
    domainD <- domainD[domainD$conScore >= 12,]
    DD <- domainD$start[2:nrow(domainD)] - domainD$end[1:(nrow(domainD)-1)]
    starts <- c(domainD[1,"start"] ,domainD[(1:nrow(domainD))[DD > 1] + 1,"start"])
    ends <- c(domainD[DD > 1,"end"], domainD[nrow(domainD),"end"])
    if(is.na(starts[length(starts)])){
      starts = starts[1:(length(starts) - 1)]
      ends = ends[1:(length(ends) - 1)]
    }
    domainRanges = rbind(domainRanges, data.frame(chr = chr, start = starts, end = ends, domain = d))
  }
}

chrName <- unique(as.character(conDomain$chr))
for(chr in chrName){
  domains <- domainRanges[domainRanges$chr == chr,]
  domains = domains[order(domains$start),]
  domainRanges[domainRanges$chr == chr,] = domains
}


allRangeE <- NULL
allRangeL <- NULL
allRange <- NULL
for(i in 1:length(sampNames)){
  ranges <- get(sampNames[i])
  ranges$sample <- sampNames[i]
  allRange <- rbind(allRange, ranges)
  allRangeE <- rbind(allRangeE, ranges[ranges$domain == "ERD",])
  allRangeL <- rbind(allRangeL, ranges[ranges$domain == "LRD",])
}
ranges <- domainRanges
ranges$sample <- "constitutive"
allRange <- rbind(allRange, ranges)
allRangeE <- rbind(allRangeE, ranges[ranges$domain == "ERD",])
allRangeL <- rbind(allRangeL, ranges[ranges$domain == "LRD",])

pdf(file ="plots/consituativeDomains/DomainSizes.pdf", width = 10, height =7)
layout(matrix(c(1,2), nrow = 2))
par(mar = c(2,5,8,5))
boxplot(log10(allRangeE$end-allRangeE$start) ~ as.factor(allRangeE$sample), las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "domain size (log10 bp)", names = rep("", 17),ylim = c(3,8), main = "Replication Domain Size")
axis(side = 4,at = 3 + ((8-3)/2), labels = "Early")
par(mar = c(10,5,0,5))
boxplot(log10(allRangeL$end-allRangeL$start) ~ as.factor(allRangeL$sample), las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "domain size (log10 bp)", ylim = c(3,8))
axis(side = 4,at = 3 + ((8-3)/2), labels = "Late")
mtext(text = "Sample", side = 1, line = 8, cex = 1.2)
dev.off()


par(mar = c(5,5,5,5))

# so on average we are getting smaller domains
# may be becasue we are cutting them when we get to an unknown bin 
domainRangesE <- domainRanges[domainRanges$domain == "ERD",]
domainRangesL <- domainRanges[domainRanges$domain == "LRD",]

# difference between boundaries
rangeTestE <- domainRangesE
differencesE <- rangeTestE$start[2:nrow(rangeTestE)] - rangeTestE$end[1:(nrow(rangeTestE)-1)]
MostRangeE <- allRangeE[allRangeE$sample != "constitutive",]
differencesALLE <- MostRangeE$start[2:nrow(MostRangeE)] - MostRangeE$end[1:(nrow(MostRangeE)-1)]

rangeTestL <- domainRangesL
differencesL <- rangeTestL$start[2:nrow(rangeTestL)] - rangeTestL$end[1:(nrow(rangeTestL)-1)]
MostRangeL <- allRangeL[allRangeL$sample != "constitutive",]
differencesALLL <- MostRangeL$start[2:nrow(MostRangeL)] - MostRangeL$end[1:(nrow(MostRangeL)-1)]


pdf(file = "plots/consituativeDomains/distanceBetweenDomain.pdf")
layout(matrix(c(1,2), nrow = 2))
plot(density(log10(differencesE[differencesE > 0])), ylim = c(0,1.2), xlim = c(2,8),
     main = "distances between adjacent early replicating domains", xlab = "distances (log10 bp)")
lines(density(log10(differencesALLE[differencesALLE > 0 ])), col = 2)
legend("topleft",legend = c("constitutive domains", "pooled domains"), fill = c(1,2))

plot(density(log10(differencesL[differencesL > 0])), ylim = c(0,1.2), xlim = c(2,8),
     main = "distances between adjacent late replicating domains", xlab = "distances (log10 bp)")
lines(density(log10(differencesALLL[differencesALLL > 0 ])), col = 2)
legend("topleft",legend = c("constitutive domains", "pooled domains"), fill = c(1,2))
dev.off()
# the gaps between adjacent ranges 


### the next step is to look at how many of our ranges are falling within domains


# ref gene
web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/refGene.txt.gz", sep = "")
con <- gzcon(url(web))
txt <- readLines(con)
refgene <- read.delim(textConnection(txt), header = FALSE)

refgene.gr <- GRanges(seqnames=Rle(refgene[,3]), ranges = IRanges(start = refgene[,5], end = refgene[,6]))
refgene_gap.gr <- filterIntergenic(refgene)
intronKeep.gr <- filterIntron(refgene)


# the proportion of intervals within 

E.gr <- GRanges(seqnames = Rle(allRangeE$chr),
                ranges = IRanges(start = allRangeE$start, end = allRangeE$end),
                domain = allRangeE$domain,
                sample = allRangeE$sample)

L.gr <- GRanges(seqnames = Rle(allRangeL$chr),
                ranges = IRanges(start = allRangeL$start, end = allRangeL$end),
                domain = allRangeL$domain,
                sample = allRangeL$sample)


####### early replicating

IntergenicOLE <- as.matrix(findOverlaps(refgene_gap.gr, E.gr, type = "within"))

IntergenicdfE <- data.frame(IntergenicBases = width(refgene_gap.gr[IntergenicOLE[,1]]), 
                  samples = as.factor(elementMetadata(E.gr)[["sample"]][IntergenicOLE[,2]]),
                  RangesBases = width(E.gr[IntergenicOLE[,2]]))

summary(IntergenicdfE$samples)
aggIntergenicE<- aggregate(x = IntergenicdfE$IntergenicBases,by = list(IntergenicdfE$samples), FUN = sum)
aggSamplesE <- aggregate(x = width(E.gr),by = list(elementMetadata(E.gr)[["sample"]]), FUN = sum)


DFmergeEintergenic <- merge(aggIntergenicE, aggSamplesE, by.x = 1, by.y = 1)
DFmergeEintergenic$portion = DFmergeEintergenic[,2]/DFmergeEintergenic[,3]
colnames(DFmergeEintergenic) <- c("sample", "intergenicBp", "ERDbp", "ERD_Intergenic_Proportion")

####Intron
IntronOLE <- as.matrix(findOverlaps(intronKeep.gr, E.gr, type = "within"))

IntrondfE <- data.frame(IntronBases = width(intronKeep.gr[IntronOLE[,1]]), 
                            samples = as.factor(elementMetadata(E.gr)[["sample"]][IntronOLE[,2]]),
                            RangesBases = width(E.gr[IntronOLE[,2]]))


aggIntronE<- aggregate(x = IntrondfE$IntronBases,by = list(IntrondfE$samples), FUN = sum)
aggSamplesE <- aggregate(x = width(E.gr),by = list(elementMetadata(E.gr)[["sample"]]), FUN = sum)


DFmergeEintron <- merge(aggIntronE, aggSamplesE, by.x = 1, by.y = 1)
DFmergeEintron$portion = DFmergeEintron[,2]/DFmergeEintron[,3]
colnames(DFmergeEintron) <- c("sample", "intronBp", "ERDbp", "ERD_Intron_Proportion")



####### late replicating
### Intergenic L

IntergenicOLL <- as.matrix(findOverlaps(refgene_gap.gr, L.gr, type = "within"))

IntergenicdfL <- data.frame(IntergenicBases = width(refgene_gap.gr[IntergenicOLL[,1]]), 
                            samples = as.factor(elementMetadata(L.gr)[["sample"]][IntergenicOLL[,2]]),
                            RangesBases = width(L.gr[IntergenicOLL[,2]]))

summary(IntergenicdfL$samples)
aggIntergenicL<- aggregate(x = IntergenicdfL$IntergenicBases,by = list(IntergenicdfL$samples), FUN = sum)
aggSamplesL <- aggregate(x = width(L.gr),by = list(elementMetadata(L.gr)[["sample"]]), FUN = sum)


DFmergeLintergenic <- merge(aggIntergenicL, aggSamplesL, by.x = 1, by.y = 1)
DFmergeLintergenic$portion = DFmergeLintergenic[,2]/DFmergeLintergenic[,3]
colnames(DFmergeLintergenic) <- c("sample", "intergenicBp", "LRDbp", "LRD_Intergenic_Proportion")


##### Intron L

IntronOLL <- as.matrix(findOverlaps(intronKeep.gr, L.gr, type = "within"))

IntrondfL <- data.frame(IntronBases = width(intronKeep.gr[IntronOLL[,1]]), 
                        samples = as.factor(elementMetadata(L.gr)[["sample"]][IntronOLL[,2]]),
                        RangesBases = width(L.gr[IntronOLL[,2]]))


aggIntronL<- aggregate(x = IntrondfL$IntronBases,by = list(IntrondfL$samples), FUN = sum)
aggSamplesL <- aggregate(x = width(L.gr),by = list(elementMetadata(L.gr)[["sample"]]), FUN = sum)


DFmergeLintron <- merge(aggIntronL, aggSamplesL, by.x = 1, by.y = 1)
DFmergeLintron$portion = DFmergeLintron[,2]/DFmergeLintron[,3]
colnames(DFmergeLintron) <- c("sample", "intronBp", "LRDbp", "LRD_Intron_Proportion")


ErdMerge <- merge(DFmergeEintergenic, DFmergeEintron, by.x = 1, by.y = 1)
LrdMerge <- merge(DFmergeLintergenic, DFmergeLintron, by.x = 1, by.y = 1)
DomainMerge <- merge(ErdMerge, LrdMerge, by.x = 1, by.y =1)

# We've got a sample of 

# basicly we want to know if we have stuff on the same order fo magnitude

# Should we put this in a table 


ordering <- boxplot(log10(IntergenicdfE$IntergenicBases) ~ IntergenicdfE$samples)
orderNo <- order(ordering$stats[3,])
orderName <- ordering$names[order(ordering$stats[3,])]

##### plotting the Data

IntergenicdfL$samples <- factor(IntergenicdfL$samples,levels = orderName)
IntergenicdfE$samples <- factor(IntergenicdfE$samples,levels = orderName)
IntrondfL$samples <- factor(IntrondfL$samples,levels = orderName)
IntrondfE$samples <- factor(IntrondfE$samples,levels = orderName)


pdf(file = "plots/consituativeDomains/IntervalSizesInDomains.pdf", onefile = T, width = 10, height = 6)


layout(matrix(c(1,2,3,4), nrow = 4),heights = c(1,1,.6,1.75))

par(mar = c(0,5,5,5))

boxplot(log10(IntrondfE$IntronBases) ~ IntrondfE$samples,
        ylab = "", names = rep("", 17), main = "Intergenic region Size",
        axes = FALSE, ylim =c(6,8))

text(1:nrow(DFmergeEintergenic),rep(7.75, nrow(DFmergeEintergenic)),labels = round(DFmergeEintergenic$intergenicBp/1000000)[orderNo])
text(1:nrow(DFmergeEintergenic),rep(7, nrow(DFmergeEintergenic)),labels = round(DFmergeEintergenic$ERDbp/1000000)[orderNo])
text(1:nrow(DFmergeEintergenic),rep(6.25, nrow(DFmergeEintergenic)),labels = round(DFmergeEintergenic$ERD_Intergenic_Proportion * 100,digits =1)[orderNo])

mtext("Intergenic\n(Mb)", side = 2,at = 7.75, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Early domain\n(Mb)", side = 2,at = 7, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Intergenic\ncoverage (%)", side = 2,at = 6.25, las = 1, line = 3.5,cex = .7,adj = 0)



par(mar = c(2,5,0,5))

boxplot(log10(IntergenicdfE$IntergenicBases) ~ IntergenicdfE$samples, las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "interval size (log10 bp)", names = rep("", 17),ylim = c(2,7), main = "",
        col = c(2, rep(0,16)))
axis(side = 4,at = 2 + ((7-2)/2), labels = "Early")
abline(h = mean(log10(IntergenicdfL$IntergenicBases)), col = 4, lty = 2)


par(mar = c(0,5,0,5))

boxplot(log10(IntrondfE$IntronBases) ~ IntrondfE$samples,
        ylab = "", names = rep("", 17),
        axes = FALSE, ylim =c(6,8))

text(1:nrow(DFmergeLintergenic),rep(7.75, nrow(DFmergeLintergenic)),labels = round(DFmergeLintergenic$intergenicBp/1000000)[orderNo])
text(1:nrow(DFmergeLintergenic),rep(7, nrow(DFmergeLintergenic)),labels = round(DFmergeLintergenic$LRDbp/1000000)[orderNo])
text(1:nrow(DFmergeLintergenic),rep(6.25, nrow(DFmergeLintergenic)),labels = round(DFmergeLintergenic$LRD_Intergenic_Proportion * 100,digits =1)[orderNo])

mtext("Intergenic\n(Mb)", side = 2,at = 7.75, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Late domain\n(Mb)", side = 2,at = 7, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Intergenic\ncoverage (%)", side = 2,at = 6.25, las = 1, line = 3.5,cex = .7,adj = 0)


par(mar = c(10,5,0,5))
boxplot(log10(IntergenicdfL$IntergenicBases) ~ IntergenicdfL$samples, las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "interval size (log10 bp)",ylim = c(2,7), main = "",col = c( 2, rep(0,16)))

axis(side = 4,at = 2 + ((7-2)/2), labels = "Late")
abline(h = mean(log10(IntergenicdfE$IntergenicBases)), col = 4, lty = 2)
mtext(text = "Sample", side = 1, line = 8, cex = 1.2)




######### Plot Intron 

layout(matrix(c(1,2,3,4), nrow = 4),heights = c(1,1,.6,1.75))

par(mar = c(0,5,5,5))

boxplot(log10(IntrondfE$IntronBases) ~ IntrondfE$samples,
        ylab = "", names = rep("", 17), main = "Intron region Size",
        axes = FALSE, ylim =c(6,8))

text(1:nrow(DFmergeEintron),rep(7.75, nrow(DFmergeEintron)),labels = round(DFmergeEintron$intronBp/1000000)[orderNo])
text(1:nrow(DFmergeEintron),rep(7, nrow(DFmergeEintron)),labels = round(DFmergeEintron$ERDbp/1000000)[orderNo])
text(1:nrow(DFmergeEintron),rep(6.25, nrow(DFmergeEintron)),labels = round(DFmergeEintron$ERD_Intron_Proportion * 100,digits =1)[orderNo])

mtext("Intron\n(Mb)", side = 2,at = 7.75, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Early domain\n(Mb)", side = 2,at = 7, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Intron coverage\n(%)", side = 2,at = 6.25, las = 1, line = 3.5,cex = .7,adj = 0)



par(mar = c(2,5,0,5))

boxplot(log10(IntrondfE$IntronBases) ~ IntrondfE$samples, las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "interval size (log10 bp)", names = rep("", 17),ylim = c(2,6), main = "",
        col = c(2, rep(0,16)))
axis(side = 4,at = 2 + ((6-2)/2), labels = "Early")
abline(h = mean(log10(IntrondfL$IntronBases)), col = 4, lty = 2)


par(mar = c(0,5,0,5))

boxplot(log10(IntrondfE$IntronBases) ~ IntrondfE$samples,
        ylab = "", names = rep("", 17),
        axes = FALSE, ylim =c(6,8))

text(1:nrow(DFmergeLintron),rep(7.75, nrow(DFmergeLintron)),labels = round(DFmergeLintron$intronBp/1000000)[orderNo])
text(1:nrow(DFmergeLintron),rep(7, nrow(DFmergeLintron)),labels = round(DFmergeLintron$LRDbp/1000000)[orderNo])
text(1:nrow(DFmergeLintron),rep(6.25, nrow(DFmergeLintron)),labels = round(DFmergeLintron$LRD_Intron_Proportion * 100,digits =1)[orderNo])

mtext("Intron\n(Mb)", side = 2,at = 7.75, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Late domain\n(Mb)", side = 2,at = 7, las = 1, line =  3.5,cex = .7,adj = 0)
mtext("Intron coverage\n(%)", side = 2,at = 6.25, las = 1, line = 3.5,cex = .7,adj = 0)


par(mar = c(10,5,0,5))
boxplot(log10(IntrondfL$IntronBases) ~ IntrondfL$samples, las = 2, varwidth = TRUE, notch = TRUE,
        ylab = "interval size (log10 bp)", ylim = c(2,6), col = c( 2, rep(0,16)))
axis(side = 4,at = 2 + ((6-2)/2), labels = "Late")
abline(h = mean(log10(IntrondfE$IntronBases)), col = 4, lty = 2)
mtext(text = "Sample", side = 1, line = 8, cex = 1.2)


dev.off()



write.table(x = rbind(domainRangesE,domainRangesL), file = "Data/ConsTimingDomains", quote = F, sep = "\t", row.names = F, col.names = T)



# so now we have a bunch of plots for the sup to support our constituative domain ID

