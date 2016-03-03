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
waveDomain[waveMean > 0 + waveSd] <- 1
waveDomain[waveMean < 0 - waveSd] <- -1
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

plot(waveMean[1:100000])
par(new=T)
plot(x = 1001:100000, y = abs(slope), type = "l", xlim = c(0,100000))


ERD <- domainRanges[domainRanges$domain == 1,]
ERD.gr <- GRanges(seqnames = Rle(ERD$chr), 
                  ranges = IRanges(start = ERD$start, end = ERD$end))

new_L1.gr <- GRanges(seqnames = Rle(joinRepChromatin$new_L1$genoName),
                     ranges = IRanges(start = joinRepChromatin$new_L1$genoStart, end = joinRepChromatin$new_L1$genoEnd))

int.ref <- intersect(ERD.gr, refgene.gr)

int <- intersect(ERD.gr,new_L1.gr)

new_L1.OL <- as.matrix(findOverlaps(ERD.gr, int))
new_L1.agg <- aggregate(x = width(int)[new_L1.OL[,2]], by = list(new_L1.OL[,1]), FUN = sum)

refgene.OL <- as.matrix(findOverlaps(ERD.gr, int.ref))
refgene.agg <- aggregate(x = width(int.ref)[refgene.OL[,2]], by = list(refgene.OL[,1]), FUN = sum)



plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log( (new_L1.agg$x)/(width(ERD.gr)[new_L1.agg$Group.1])))
plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log(new_L1.agg$x))
abline(lm(log(new_L1.agg$x) ~ log(width(ERD.gr)[new_L1.agg$Group.1])))

## so far the smaller the bin the higher the L1 density

plot(log(width(ERD.gr)[new_L1.agg$Group.1]),log( (new_L1.agg$x)/(width(ERD.gr)[new_L1.agg$Group.1])))

plot(log(width(ERD.gr)[refgene.agg$Group.1]),log(refgene.agg$x))



A <- merge(new_L1.agg, refgene.agg, by = 1)
plot(    ( (A$x.y)/(width(ERD.gr)[A$Group.1]))   ,   ( (A$x.x)/(width(ERD.gr)[A$Group.1])))
cor(    ( (A$x.y)/(width(ERD.gr)[A$Group.1]))   ,   ( (A$x.x)/(width(ERD.gr)[A$Group.1])))


plot( log10( (A$x.y))   ,  log10 ( (A$x.x)))



