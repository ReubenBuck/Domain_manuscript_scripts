# lets get the gene expression info working,

# quite simple, we stratify our genes by expresssion level in each sample 
# then at each expression level we plot the inton size distribution, 
# the upstream size distribution 
# the downstream size distribution 




setwd("~/Desktop/Domain_manuscript/")

rm(list = ls())


library(GenomicRanges)
library(rtracklayer)
library(zoo)

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





tissue <- c("adipose", "brain", "liver", "skeletal_muscle", "testes", "kidney")
for(i in tissue){
  expr <- read.table(paste("~/Desktop/expresionLu/",i,"_expression_output", sep =""))
  expr <- data.frame(ID = expr[,4], chr = expr[,1], start = expr[,2], stop = expr[,3], length = expr[,5], effective_length = expr[,6], TPM = expr[,8], FPKM=expr[,9])
  dup <- as.character(unique(expr$ID[duplicated(expr$ID)]))
  expr <- expr[!(expr$ID %in% dup),]
  expr1 <- expr[,1:6]
  colnames(expr)[7] <- paste(i, "_TPM", sep = "")
  colnames(expr)[8] <- paste(i, "_FPKM", sep = "")
  expr <- expr[,c(1,7,8)] 
  expr[,1] <- as.character(expr[,1])
  expr[,2] <- log2(expr[,2] )
  expr[,3] <- log2(expr[,3] )
  assign(i, expr)
}


for(i in tissue){
  expr1 <- merge(expr1,get(i), by.x=1, by.y=1)
}


Expr <- expr1
Expr <- Expr[rowSums(is.infinite(as.matrix(Expr[,6:ncol(Expr)]))) < 4,]
Expr$mean <- apply(X = Expr[,6:12],MARGIN = 1,FUN = median)
Expr$groups <- cut((Expr$mean), breaks = 10)

boxplot(log10(Expr$stop - Expr$start) ~ Expr$groups, notch = F,outline = T)

genes <- data.frame(chr = Expr$chr, start = Expr$start, end = Expr$stop)
genes.gr <- GRanges(seqnames = Rle(genes$chr), 
                    ranges = IRanges(start = genes$start, end = genes$end))
f.OL <- as.matrix(findOverlaps(genes.gr, intronKeep.gr))
ILengths <- data.frame(intron_length = width(intronKeep.gr)[f.OL[,2]], groups = Expr$groups[f.OL[,1]])

layout(matrix(c(1,2), nrow = 1))
par(mar = c(10,5,5,2))

boxplot(log10(ILengths$intron_length) ~ ILengths$groups  , outline = T, notch = T, las = 2, main = "Intronic",  ylab = "log10 length (bp)")
# very long introns tend to belong to genes with lower expression 

f.OL <- as.matrix(findOverlaps(genes.gr, refgene_gap.gr, maxgap = 3))
GapLengths <- data.frame(intergenic_length = width(refgene_gap.gr)[f.OL[,2]], groups = Expr$groups[f.OL[,1]])

par(mar = c(10,2,5,5))
boxplot(log(GapLengths$intergenic_length) ~GapLengths$groups, outline = T, notch = T, las = 2,  main = "Intergenic",  ylab = "log10 length (bp)")






