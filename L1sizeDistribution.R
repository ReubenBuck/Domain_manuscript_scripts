### so what do we need?

### reps

### intergenic regions and intronic regions 

### Domains 


# so we want to see if the distribtuion of L1 sizes changes in different regions, 
# we want to see the element size argument

# we also want to calculate which portion of the elements are above and below certain sizes. 




setwd("~/Desktop/Domain_manuscript/")

rm(list = ls())

library(GenomicRanges)
library(rtracklayer)

spec1 <- "Human"
genome = "hg19"


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/rep_db.R")

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


# get repeat info
rep = rep_info(spec1=spec1,genome=genome)
# sort into intergenic bins and intronic bins
intergenic_reps <- binSort(rep=rep, bins=bins_gene_gap, TE.names=names(rep), repType = rep("repeats",length(rep)))
intron_reps <- binSort(rep=rep, bins=bins_intron, TE.names=names(rep), repType = rep("repeats",length(rep)))




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



hist(joinRep$new_L1$repEnd, breaks = 1000)



new_L1repstart <- joinRep$new_L1$repStart
new_L1repstart[joinRep$new_L1$strand == "-"] <- joinRep$new_L1$repLeft[joinRep$new_L1$strand == "-"]
new_L1repstart[joinRep$new_L1$repEnd - new_L1repstart < 0] <- joinRep$new_L1$repEnd[joinRep$new_L1$repEnd - new_L1repstart < 0]

new_L1repleft <- joinRep$new_L1$repLeft
new_L1repleft[joinRep$new_L1$strand == "-"] <- joinRep$new_L1$repStart[joinRep$new_L1$strand == "-"]
new_L1repleft <- abs(new_L1repleft)

new_L1.gr <- GRanges(seqnames = joinRep$new_L1$genoName, 
                     ranges = IRanges(start = joinRep$new_L1$genoStart, end = joinRep$new_L1$genoEnd),
                     repRange = IRanges(end = (new_L1repleft + joinRep$new_L1$repEnd) - new_L1repstart, start = (new_L1repleft + joinRep$new_L1$repEnd) - joinRep$new_L1$repEnd))


#new_L1.gr <- new_L1.gr[grep(pattern = "L1PA",x = joinRep$new_L1$repName) ]


pdf(file = "plots/new_L1_length_analysis.pdf", onefile = TRUE)

new_L1cov <- as.integer(coverage(elementMetadata(new_L1.gr)[["repRange"]]))
plot(new_L1cov, type = "l", main = "new_L1 coverage", xlab = "position from 3' end (bp)", ylab = "coverage")
lo <- predict(loess(formula = new_L1cov ~ seq(1,length(new_L1cov)),span = .1),se = F)
lines(lo, col = 2)
#lines(lo$fit+3*lo$se.fit, col = 2)
#lines(lo$fit-3*lo$se.fit, col = 2)
legend("topright", c("observed", "model fit"), lty=1, col = c(1,2))

new_L1cov <- lo


## we need our four groups 
## L1s intergenic E
## L1s intergenic L
## L1s intronic E
## L1s intornic L

new_L1intergenicE <- subsetByOverlaps(subsetByOverlaps(new_L1.gr, domainE.gr, type = "within"),refgene_gap.gr, type = "within")
new_L1intergenicEcov <- as.integer(coverage(elementMetadata(new_L1intergenicE)[[1]]))
new_L1intergenicElo <- predict(loess(formula = new_L1intergenicEcov ~ seq(1,length(new_L1intergenicEcov)),span = .1),newdata = seq(1,length(new_L1intergenicEcov)))

new_L1intergenicL <- subsetByOverlaps(subsetByOverlaps(new_L1.gr, domainL.gr, type = "within"),refgene_gap.gr, type = "within")
new_L1intergenicLcov <- as.integer(coverage(elementMetadata(new_L1intergenicL)[[1]]))
new_L1intergenicLlo <- predict(loess(formula = new_L1intergenicLcov ~ seq(1,length(new_L1intergenicLcov)),span = .1),newdata = seq(1,length(new_L1intergenicLcov)))

new_L1intronE <- subsetByOverlaps(subsetByOverlaps(new_L1.gr, domainE.gr, type = "within"),intronKeep.gr, type = "within")
new_L1intronEcov <- as.integer(coverage(elementMetadata(new_L1intronE)[[1]]))
new_L1intronElo <- predict(loess(formula = new_L1intronEcov ~ seq(1,length(new_L1intronEcov)),span = .1),newdata = seq(1,length(new_L1intronEcov)))

new_L1intronL <- subsetByOverlaps(subsetByOverlaps(new_L1.gr, domainL.gr, type = "within"),intronKeep.gr, type = "within")
new_L1intronLcov <- as.integer(coverage(elementMetadata(new_L1intronL)[[1]]))
new_L1intronLlo <- predict(loess(formula = new_L1intronLcov ~ seq(1,length(new_L1intronLcov)),span = .1),newdata = seq(1,length(new_L1intronLcov)))



plot(new_L1intergenicElo / sum(width(new_L1intergenicE)), type = "l", lwd = 2, ylim = c(0,.001),
     xlab = "position form L1 3' end", ylab = "coverage density", main = "element length estimation", col = 2)
lines(new_L1intergenicLlo / sum(width(new_L1intergenicL)), col = 3, lwd = 2)

lines(new_L1intronElo / sum(width(new_L1intronE)), type = "l", lty = 3, lwd = 2, col =2)
lines(new_L1intronLlo / sum(width(new_L1intronL)), col = 3, lty = 3, lwd = 2)

legend("topright", legend = c("intergenic early", "intergenic late", "intron early", "intron late"), lty = c(1,1,3,3), col = c(2,3,2,3), lwd = c(2,2,2,2))


ratioIntergenic <- as.double(new_L1intergenicElo / sum(width(new_L1intergenicE)))[1:6000]/
  as.double(new_L1intergenicLlo / sum(width(new_L1intergenicL)))[1:6000]
  
ratioIntron <- as.double(new_L1intronElo / sum(width(new_L1intronE)))[1:6000]/
  as.double(new_L1intronLlo / sum(width(new_L1intronL)))[1:6000]


  
plot(ratioIntergenic, type = "l", ylim = c(.5,1.5), lwd =2, ylab = "Early:Late coverage density", xlab = "position form L1 3' end", 
     main = "estimated element length ratio")  
lines(ratioIntron, lty = 3, lwd = 2)
legend("topright", legend = c("intergenic", "intron"), lty = c(1,3), lwd = c(2,2))
abline(h = 1, col = 4)

# so now we have a plot of ratios

# next setp is to look at how may elememts there are at a certain size 

# maybe we set up bins of 100 bp? 
# we can do the cummulative density and just minus the part of the curve we are not interested in

cdfRate <- function(len, cov){
  sum(cov[1:len] - cov[len])/
    sum(cov)
}

cdfSum <- function(len, cov){
  sum(cov[1:len] - cov[len])
}



newL1intronEcD <- sapply(X = 1:length(new_L1intronElo), FUN = cdfSum, cov = new_L1intronElo)
newL1intronLcD <- sapply(X = 1:length(new_L1intronLlo), FUN = cdfSum, cov = new_L1intronLlo)

newL1intergenicEcD <- sapply(X = 1:length(new_L1intergenicElo), FUN = cdfSum, cov = new_L1intergenicElo)
newL1intergenicLcD <- sapply(X = 1:length(new_L1intergenicLlo), FUN = cdfSum, cov = new_L1intergenicLlo)


intronEbases <- sum(width(subsetByOverlaps(intronKeep.gr,domainE.gr, type = "within")))
intronLbases <- sum(width(subsetByOverlaps(intronKeep.gr,domainL.gr, type = "within")))
intergenicEbases <- sum(width(subsetByOverlaps(refgene_gap.gr,domainE.gr, type = "within")))
intergenicLbases <- sum(width(subsetByOverlaps(refgene_gap.gr,domainL.gr, type = "within")))

sizeRange <- 50

plot((newL1intergenicLcD[seq(sizeRange +1, length(newL1intergenicLcD), sizeRange)] - 
        newL1intergenicLcD[seq(1, length(newL1intergenicLcD) - sizeRange, sizeRange)])/intergenicLbases , 
     xaxt = "n", type = "l", col = 3, lwd = 2, xlab = "length from 3' end (bp)",
     ylab = "region coverage per pb", main = paste("L1s between length x and x +",sizeRange, "bp"))
lines((newL1intergenicEcD[seq(sizeRange +1, length(newL1intergenicEcD), sizeRange)] - 
         newL1intergenicEcD[seq(1, length(newL1intergenicEcD) - sizeRange, sizeRange)])/intergenicEbases , 
      xaxt = "n", type = "l", col = 2, lwd = 2)

lines((newL1intronLcD[seq(sizeRange +1, length(newL1intronLcD), sizeRange)] - 
         newL1intronLcD[seq(1, length(newL1intronLcD) - sizeRange, sizeRange)])/intronLbases , 
      xaxt = "n", type = "l", col = 3, lty = 3, lwd = 2)
lines((newL1intronEcD[seq(sizeRange +1, length(newL1intronEcD), sizeRange)] - 
         newL1intronEcD[seq(1, length(newL1intronEcD) - sizeRange, sizeRange)])/intronEbases , 
      xaxt = "n", type = "l", col = 2, lty = 3, lwd = 2)

axis(side = 1,at = seq(0,length(newL1intronLcD), by = 10),labels = seq(0,length(newL1intronLcD), by = 10)*sizeRange)
legend("topright", legend = c("intergenic early", "intergenic late", "intron early", "intron late"), lty = c(1,1,3,3), col = c(2,3,2,3), lwd = c(2,2,2,2))


intergenicRatioRanges <- (((newL1intergenicEcD[seq(sizeRange +1, length(newL1intergenicEcD), sizeRange)] - 
                              newL1intergenicEcD[seq(1, length(newL1intergenicEcD) - sizeRange, sizeRange)])/intergenicEbases )/
                            ((newL1intergenicLcD[seq(sizeRange +1, length(newL1intergenicLcD), sizeRange)] - 
                             newL1intergenicLcD[seq(1, length(newL1intergenicLcD) - sizeRange, sizeRange)])/intergenicLbases))
intronRatioRanges <- (((newL1intronEcD[seq(sizeRange +1, length(newL1intronEcD), sizeRange)] - 
                          newL1intronEcD[seq(1, length(newL1intronEcD) - sizeRange, sizeRange)])/intronEbases )/
                        ((newL1intronLcD[seq(sizeRange +1, length(newL1intronLcD), sizeRange)] - 
                            newL1intronLcD[seq(1, length(newL1intronLcD) - sizeRange, sizeRange)])/intronLbases))

plot(intergenicRatioRanges, type = "l", lwd = 2, xaxt = "n", ylim = c(0,2), 
     xlab = "length from 3' end (bp)",
     ylab = "E:L region coverage per bp", main = paste("L1s between length x and x +",sizeRange, "bp"))
lines(intronRatioRanges, lty = 3, lwd = 2)
axis(side = 1,at = seq(0,length(newL1intronLcD), by = 10),labels = seq(0,length(newL1intronLcD), by = 10)*sizeRange)
legend("topright", legend = c("intergenic", "intron"), lty = c(1,3), lwd = c(2,2))




plot(newL1intergenicLcD /intergenicLbases , 
      type = "l", col = 3, lwd = 2, xlab = "length from 3' end (bp)",
     ylab = "region coverage per pb", main = "all L1s < x bp")
lines(newL1intergenicEcD/intergenicEbases , 
       type = "l", col = 2, lwd = 2)

lines(newL1intronLcD/intronLbases , 
       type = "l", col = 3, lty = 3, lwd = 2)
lines(newL1intronEcD/intronEbases , 
       type = "l", col = 2, lty = 3, lwd = 2)
legend("topleft", legend = c("intergenic early", "intergenic late", "intron early", "intron late"), lty = c(1,1,3,3), col = c(2,3,2,3), lwd = c(2,2,2,2))


plot((newL1intergenicEcD/intergenicEbases)[1:8000]/(newL1intergenicLcD /intergenicLbases)[1:8000], 
     type = "l", col = 1, lwd = 2, ylim = c(0,2), xlab = "length from 3' end (bp)",
     ylab = "E:L region coverage per pb", main = "all L1s < x bp")
lines((newL1intronEcD/intronEbases)/(newL1intronLcD/intronLbases) , 
      type = "l", col = 1, lty = 3, lwd = 2)
legend("topright", legend = c("intergenic", "intron"), lty = c(1,3), lwd = c(2,2))

abline(h = 1, col = 2)

dev.off()
# the conclusion from this is that size of L1s has only a small effect on their distribution
# the difference in distribution might be related to some other quality. 


head(joinRep$Alu)
AluRange <- GRanges(seqnames = Rle(joinRep$Alu$genoName), 
                    ranges = IRanges(start = joinRep$Alu$genoStart, end = joinRep$Alu$genoEnd),
                    repCoor = IRanges(start = joinRep$Alu$repStart, end = joinRep$Alu$repEnd)
)

AluE <- subsetByOverlaps(AluRange, domainE.gr)
AluL <- subsetByOverlaps(AluRange, domainL.gr)
AluEcov <- coverage(elementMetadata(AluE)$repCoor)
AluLcov <- coverage(elementMetadata(AluL)$repCoor)

plot(AluEcov,type = "l" , ylim = c(0, 600000))
lines(AluLcov,col = 2)
abline(h = length(AluE))
abline(h = length(AluL))

max(AluEcov)/length(AluE)
max(AluLcov)/length(AluL)

length(AluE)/max(AluEcov)
length(AluL)/max(AluLcov)

plot(AluEcov/sum(width(AluE)), type = "l")
lines(AluLcov/sum(width(AluL)), col= 2)

plot( (AluEcov/sum(width(AluE))) / (AluLcov/sum(width(AluL))), type = "l")

# we have the potential here to create new 


# so maybe if we join the two regions intergenic and intronic 

plot((new_L1intergenicLlo + new_L1intronElo)/((sum(width(new_L1intergenicE))) + sum(width(new_L1intronE))))
plot((new_L1cov[1:6000]), type = "l" , ylab = "", xlab = "", yaxt = "n", ylim = c(0,max(new_L1cov)), lwd = 3)



ratio <- as.double((new_L1intergenicElo[1:6000] + new_L1intronElo[1:6000]) / (sum(width(new_L1intergenicE)) + sum(width(new_L1intronE))))/
  as.double((new_L1intergenicLlo[1:6000] + new_L1intronLlo[1:6000]) / (sum(width(new_L1intergenicL)) + sum(width(new_L1intronL))))
  
pdf(file = "plots/TEopenChromInteract/L1cov.pdf", width = 3, height = 7)
layout(c(1,2))

#plot((new_L1cov[1:6000]), type = "l" , ylab = "new_L1 coverage (x1000)", xlab = "", yaxt = "n", ylim = c(0,max(new_L1cov)), lwd = 3)
plot((new_L1cov[1:6000]), type = "l" , ylab = "", xlab = "", yaxt = "n", ylim = c(0,max(new_L1cov)), lwd = 3)
axis(2,at = seq(0,max(new_L1cov), by = 25000), las = 2, labels = seq(0,max(new_L1cov), by = 25000)/1000)


plot((ratio), type = "l", ylab = "", xlab = "", ylim = c(.6,1.4), lwd = 3)
#plot((ratio), type = "l", ylab = "Early:Late coverage density", xlab = "element position (5' - 3')", ylim = c(.6,1.4), lwd = 3)
#abline(h = 1, , lty = 2, lwd = 3)
dev.off()

pdf(file = "plots/TEopenChromInteract/legends.pdf")
plot.new()
legend("center", legend = c("cERD", "cLRD"), fill = c("red", "green"), title = "Domain",box.col = "white")
legend("bottom", legend = c("exon", "non-exon"), lty = c(3,1), lwd = 3, title = "DNase1 cluster", box.col = "white")
dev.off()
