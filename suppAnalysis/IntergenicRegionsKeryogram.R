#### Intergenic region Locations

### We just need to get a keryogram 

rm(list=ls())

setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")


genome = "hg19"
spec1 = "Human"


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

plot(density(log10(bins_gene_gap$end - bins_gene_gap$start)))

hist(log10(bins_gene_gap$end - bins_gene_gap$start), breaks = 1000, xlim = c(2,7))
hist(log10(bins_intron$end - bins_intron$start), breaks = 1000, xlim = c(2,7))



chromNames <- gsub(pattern = "chr", replacement = "",x = chrom_info$V1)
chromLengths <- chrom_info[,2]
names(chromLengths) <- chromNames
chromLengths <- chromLengths[c(as.character(1:22),"X", "Y")]



longInterGenic <- bins_gene_gap[order(bins_gene_gap$end - bins_gene_gap$start,decreasing = TRUE),][1:300,]
longInterGenic <- bins_gene_gap[bins_gene_gap$end - bins_gene_gap$start > 1000000,]

longInterGenic$chr <- gsub(pattern = "chr", replacement = "", x = longInterGenic$chr)
longInterGenic$chr[longInterGenic$chr == "X"] <- 23
longInterGenic$chr[longInterGenic$chr == "Y"] <- 24
longInterGenic$chr <- as.integer(longInterGenic$chr)


longGaps <- gaps[gaps$type == "centromere" | gaps$type == "short_arm" | gaps$type == "telomere" |  gaps$type =="heterochromatin",]
longGaps$chrom <- gsub(pattern = "chr", replacement = "", longGaps$chrom)
longGaps$chrom[longGaps$chrom == "X"] <- 23
longGaps$chrom[longGaps$chrom == "Y"] <- 24
longGaps$chrom <- as.integer(longGaps$chrom)
longGaps <- longGaps[!is.na(longGaps$chrom),]



pdf(file = "plots/consituativeDomains/LongIntergenicRegions.pdf", onefile = T, width = 10, height = 6)
plot(c(1,length(chromLengths)), c(1,max(chromLengths)), frame.plot = FALSE, 
     axes = F, xlab = "chromosome", ylab = "Mb", type = "n", main = "Long intergenic regions")
axis(side = 1, at = c(1:length(chromLengths)), labels = names(chromLengths))
axis(side = 2, at = seq(0,max(chromLengths)+50000000,50000000), labels = seq(0,max(chromLengths)+50000000,50000000) / 1000000)

for(i in 1:length(chromLengths)){
  rect(xleft = i - .35, ybottom = 0, xright = i +.35, ytop = chromLengths[i])
}

for(i in 1:nrow(longInterGenic)){
  rect(xleft = longInterGenic[i,"chr"] - .35, 
       ybottom = longInterGenic[i,"start"], 
       xright = longInterGenic[i,"chr"], 
       ytop = longInterGenic[i,"end"], col = 2,
       lty = 0
       )
}



for(i in 1:nrow(longGaps)){
  rect(xleft = longGaps[i,"chrom"], 
       ybottom = longGaps[i,"chromStart"], 
       xright = longGaps[i,"chrom"] +.35, 
       ytop = longGaps[i,"chromEnd"], col = as.integer(longGaps[i,"type"]),
       lty = 0
  )
}

legend(16,max(chromLengths),legend = c("interval > 1 Mb",unique(as.character(longGaps$type))), fill = c(2,4,1,6,5))

dev.off()
