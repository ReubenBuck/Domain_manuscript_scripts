# Here we put the final nail in the coffin interms of figures 







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

Intron.GR <- GRanges(seqnames = Rle(intron_reps$counts$chr),
                     ranges = IRanges(start = intron_reps$counts$start,
                                      end = intron_reps$counts$end)
)

OLintergenicE <- as.matrix(findOverlaps(Intergenic.GR, domainE.gr, type = "within"))
IntergenicE <- Intergenic.GR[OLintergenicE[,1]]

OLintergenicL <- as.matrix(findOverlaps(Intergenic.GR, domainL.gr, type = "within"))
IntergenicL <- Intergenic.GR[OLintergenicL[,1]]


OLintronE <- as.matrix(findOverlaps(Intron.GR, domainE.gr, type = "within"))
IntronE <- Intron.GR[OLintronE[,1]]

OLintronL <- as.matrix(findOverlaps(Intron.GR, domainL.gr, type = "within"))
IntronL <- Intron.GR[OLintronL[,1]]




# final touch to the paper in figure form 
# we plot genes in each region 
# intergenic and intron
# early and late

EintergenicCounts <- intergenicJoinedFam[OLintergenicE[,1],]
LintergenicCounts <- intergenicJoinedFam[OLintergenicL[,1],]

EintronCounts <- intronJoinedFam[OLintronE[,1],]
LintronCounts <- intronJoinedFam[OLintronL[,1],]



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



#OpenCell <- allOpen$H1_ES
#OpenCell <- OpenCell[OpenCell$ocType == "DNaseOnly",]
nonExonRegion.gr <- GRanges(seqnames = Rle(c(as.character(bins_gene_gap$chr),as.character(bins_intron$chr))),
                            ranges = IRanges(start = c(bins_gene_gap$start,bins_intron$start),
                                             end = c(bins_gene_gap$end,bins_intron$end)))


for(i in c(1:2)){
  
  OpenCell <- allCell
  OpenCell <- OpenCell[OpenCell$dnasePval > 1.3,]
  Open.gr <- GRanges(seqnames = Rle(OpenCell$chr),
                     ranges = IRanges(start = OpenCell$start, end = OpenCell$end))
  if(i == 1){
    Open.gr <- subsetByOverlaps(Open.gr,nonExonRegion.gr,type = "within")
  }
  if(i == 2){
    Open.gr <- subsetByOverlaps(Open.gr,gaps(nonExonRegion.gr))
  }
  
  Open.gr <- reduce(Open.gr)
  
  OpenChromList = list(openChrom = data.frame(genoName = seqnames(Open.gr), genoStart = start(Open.gr), genoEnd = end(Open.gr)))
  
  
  #OpenIntergenic.gr <- subsetByOverlaps(Open.gr,Intergenic.GR, type = "within")
  OpenIntergenic.gr <- intersect(Open.gr,Intergenic.GR)
  openIntergenicOL <- as.matrix(findOverlaps(Intergenic.GR, OpenIntergenic.gr))
  intergenicOpen <- intergenicJoinedFam[,c("chr", "start", "end","Known")]
  
  openIntergenicAgg <- aggregate(x = width(OpenIntergenic.gr[openIntergenicOL[,2]]), by = list(openIntergenicOL[,1]), FUN = sum)
  intergenicOpen[openIntergenicAgg$Group.1,"openChrom"] <- openIntergenicAgg$x
  intergenicOpen$openChrom[is.na(intergenicOpen$openChrom)] <- 0
  
  
  #OpenIntron.gr <- subsetByOverlaps(Open.gr,Intron.GR, type = "within")
  OpenIntron.gr <- intersect(Open.gr,Intron.GR)
  openIntronOL <- as.matrix(findOverlaps(Intron.GR, OpenIntron.gr))
  intronOpen <- intronJoinedFam[,c("chr", "start", "end","Known")]
  
  openIntronAgg <- aggregate(x = width(OpenIntron.gr[openIntronOL[,2]]), by = list(openIntronOL[,1]), FUN = sum)
  intronOpen[openIntronAgg$Group.1,"openChrom"] <- openIntronAgg$x
  intronOpen$openChrom[is.na(intronOpen$openChrom)] <- 0
  if(i == 1){
    NexonintronOpen <- intronOpen
    NexonintergenicOpen <- intergenicOpen
    NexonOpenChromList <- OpenChromList
  }
  if(i == 2){
    PromintronOpen <- intronOpen
    PromintergenicOpen <- intergenicOpen
    PromOpenChromList <- OpenChromList
  }
  
}


ylims = c(.15,.2,.3,.15)

region = "intron"
#TEfam = "new_L1"
pdf(file = paste("plots/TEopenChromInteract/geneRepTime", region, ".pdf", sep=""), height = 12,width = 6)
 layout(matrix(1:10, nrow = 5, byrow = TRUE))
  par(mar=c(5,5,5,5))
for(i in 1:4){
  TEfam = names(joinRep)[i]
  
  
  posStatBins = get(paste("E",region,"Counts",sep = ""))
  lenChoice = max(posStatBins$end - posStatBins$start)/2
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")
  
  TEs_posStatsE <- TEs_posStats
  
  posStatBins = get(paste("L",region,"Counts",sep = ""))
  lenChoice = max(posStatBins$end - posStatBins$start)/2
  
  TEs_posStats <- covCalcPlot5prime3prime(lenChoice = lenChoice,repChoice = TEfam,repBins = posStatBins,repList = joinRep,refgene = refgene,type = region,repType = "repeats")
  
  TEs_posStatsL <- TEs_posStats
  

 
  len <- 100000

  
  prime5E <- rev(TEs_posStatsE$rawRepCov5/TEs_posStatsE$baseFreq5prime)[1:len]
  prime5L <- rev(TEs_posStatsL$rawRepCov5/TEs_posStatsL$baseFreq5prime)[1:len]
  
  prime3E <- (TEs_posStatsE$rawRepCov3/TEs_posStatsE$baseFreq3prime)[1:len]
  prime3L <- (TEs_posStatsL$rawRepCov3/TEs_posStatsL$baseFreq3prime)[1:len]
  
  breaker <- seq(0,log10(len)+1,by = .08)
  cuter <- cut(log10(1:len),breaks = breaker )
  agg5E <- aggregate(prime5E, list(cuter), mean)
  agg5L <- aggregate(prime5L, list(cuter), mean)
  breakDF5 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
  breakDF5 <- merge(breakDF5,agg5E, by.x = 1, by.y = 1)
  breakDF5 <- merge(breakDF5,agg5L, by.x = 1, by.y = 1)
  
  agg3E <- aggregate(prime3E, list(cuter), mean)
  agg3L <- aggregate(prime3L, list(cuter), mean)
  breakDF3 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
  breakDF3 <- merge(breakDF3,agg3E, by.x = 1, by.y = 1)
  breakDF3 <- merge(breakDF3,agg3L, by.x = 1, by.y = 1)
  
  plot((breakDF5[,2]),(breakDF5[,3]),type = "l", col = 2, xlim = c(log10(len),0), lwd = 3, xlab = region, ylab = TEfam, ylim = c(0,ylims[i]), las = 2, xaxt = "n")
  lines((breakDF5[,2]),breakDF5[,4],type = "l", col = 3, lwd = 3)
  grid()
  
  if(i == 4){
    axis(side = 1, at = 0:8, 10^(0:8)/1000)
  }
  
  plot((breakDF3[,2]),breakDF3[,3],type = "l", col = 2, xlim = c(0, log10(len)), lwd = 3, xlab = region, ylab = TEfam, ylim = c(0,ylims[i]), las = 2, xaxt = "n", yaxt = "n")
  lines((breakDF3[,2]),breakDF3[,4],type = "l", col = 3, lwd = 3)
  grid()
  
  if(i == 4){
    axis(side = 1, at = 0:8, 10^(0:8)/1000)
  }
  
  
}



### plotting the open chromatin



openStyle <- c("Nexon", "Prom")
#regionType <- c("intergenic", "intron")

#for(i in 1:2){
 # region = regionType[i]
  
  for(j in 1:2){
    openS = openStyle[j]
    OpenStatBins = get(paste(openS,region,"Open", sep = ""))
    OpenChromList <- get(paste(openS,"OpenChromList", sep = ""))
    
    
    OpenStatBinsE = OpenStatBins[get(paste("OL",region,"E", sep = ""))[,1],]
    lenChoice = 100000
    Open_posStatE = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "openChrom", 
                                        repBins = OpenStatBinsE, repList = OpenChromList, 
                                        refgene = refgene, type = region, repType = "repeats")
    
    OpenStatBinsL = OpenStatBins[get(paste("OL",region,"L", sep = ""))[,1],]
    Open_posStatL = covCalcPlot5prime3prime(lenChoice = lenChoice, repChoice = "openChrom", 
                                            repBins = OpenStatBinsL, repList = OpenChromList, 
                                            refgene = refgene, type = region, repType = "repeats")
   
    assign(paste("Open_posStatE",openS, sep=""),Open_posStatE)
    assign(paste("Open_posStatL",openS, sep=""),Open_posStatL) 
    
  }
    
    # draw lines here

    
    prime5E <- rev(Open_posStatENexon$rawRepCov5/Open_posStatE$baseFreq5prime)[1:len]
    prime5L <- rev(Open_posStatLNexon$rawRepCov5/Open_posStatL$baseFreq5prime)[1:len]
    
    breaker <- seq(0,log10(len)+1,by = .08)
    cuter <- cut(log10(1:len),breaks = breaker )
    agg5E <- aggregate(prime5E, list(cuter), mean)
    agg5L <- aggregate(prime5L, list(cuter), mean)
    breakDF5 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
    breakDF5 <- merge(breakDF5,agg5E, by.x = 1, by.y = 1)
    breakDF5 <- merge(breakDF5,agg5L, by.x = 1, by.y = 1)
    
    
    plot((breakDF5[,2]),(breakDF5[,3]),type = "l", col = 2, xlim = c(log10(len),0), lwd = 3, xlab = TEfam, ylab = region, ylim = c(0,.9), lty = 1, xaxt = "n", las = 2)
    lines((breakDF5[,2]),breakDF5[,4],type = "l", col = 3, lwd = 3, lty = 1)
    grid()
    
    prime5E <- rev(Open_posStatEProm$rawRepCov5/Open_posStatEProm$baseFreq5prime)[1:len]
    prime5L <- rev(Open_posStatLProm$rawRepCov5/Open_posStatLProm$baseFreq5prime)[1:len]
    
    breaker <- seq(0,log10(len)+1,by = .08)
    cuter <- cut(log10(1:len),breaks = breaker )
    agg5E <- aggregate(prime5E, list(cuter), mean)
    agg5L <- aggregate(prime5L, list(cuter), mean)
    breakDF5 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
    breakDF5 <- merge(breakDF5,agg5E, by.x = 1, by.y = 1)
    breakDF5 <- merge(breakDF5,agg5L, by.x = 1, by.y = 1)
    
    lines((breakDF5[,2]),(breakDF5[,3]),type = "l", col = 2, lwd = 3, lty = 3)
    lines((breakDF5[,2]),breakDF5[,4],type = "l", col = 3, lwd = 3, lty = 3)
    grid()
    
    
    
  ##### 3 prime side
    
    prime3E <- (Open_posStatENexon$rawRepCov3/Open_posStatENexon$baseFreq3prime)[1:len]
    prime3L <- (Open_posStatLNexon$rawRepCov3/Open_posStatLNexon$baseFreq3prime)[1:len]
    
    breaker <- seq(0,log10(len)+1,by = .08)
    cuter <- cut(log10(1:len),breaks = breaker )

    agg3E <- aggregate(prime3E, list(cuter), mean)
    agg3L <- aggregate(prime3L, list(cuter), mean)
    breakDF3 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
    breakDF3 <- merge(breakDF3,agg3E, by.x = 1, by.y = 1)
    breakDF3 <- merge(breakDF3,agg3L, by.x = 1, by.y = 1)
    
    
    plot((breakDF3[,2]),breakDF3[,3],type = "l", col = 2, xlim = c(0, log10(len)), lwd = 3, xlab = TEfam, ylab = region, ylim = c(0,.9), lty = 1, xaxt = "n", yaxt = "n")
    lines((breakDF3[,2]),breakDF3[,4],type = "l", col = 3, lwd = 3, lty = 1)
    grid()
    
    prime3E <- (Open_posStatEProm$rawRepCov3/Open_posStatEProm$baseFreq3prime)[1:len]
    prime3L <- (Open_posStatLProm$rawRepCov3/Open_posStatLProm$baseFreq3prime)[1:len]
    
    breaker <- seq(0,log10(len)+1,by = .08)
    cuter <- cut(log10(1:len),breaks = breaker )
    
    agg3E <- aggregate(prime3E, list(cuter), mean)
    agg3L <- aggregate(prime3L, list(cuter), mean)
    breakDF3 <- data.frame(level = levels(cuter), mids = breaker[1:(length(breaker) - 1)] + .015)
    breakDF3 <- merge(breakDF3,agg3E, by.x = 1, by.y = 1)
    breakDF3 <- merge(breakDF3,agg3L, by.x = 1, by.y = 1)
    
    
    lines((breakDF3[,2]),breakDF3[,3],type = "l", col = 2, lwd = 3, lty = 3)
    lines((breakDF3[,2]),breakDF3[,4],type = "l", col = 3, lwd = 3, lty = 3)
    grid()
    
    
dev.off()
#}










