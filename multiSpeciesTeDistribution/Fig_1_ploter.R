## lets do PCA analysis 



library(ggplot2)
library(ks)
library(gplots)
library(GenomicRanges)
require(ggbio)


setwd("~/Desktop/Dist_matrix_TE_div/")

source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")


spec1 <- "Human"
spec2 <- "Horse"
spec3 <- "Bovine"
spec4 <- "Dog"
spec5 <- "Mouse"



ucsc1 <- "hg19"
rem.un <- "yes"
no.xy <- F

keep.limit <- 1350000
mb <- 1000000
# bin.size is set so that the indexing stays consitent across pipeline
bin.size = 500000

# trying out the unscaled method
keep.NGenes = "yes"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "no"

#create objects into whcih i will store the binsizes

# make a loop here to get all the species' names 


s1name <- paste("count_tables/",spec1, "_AllBinCounts.txt", sep = "")
s2name <- paste("count_tables/",spec2, "_AllBinCounts.txt", sep = "")
s3name <- paste("count_tables/",spec3, "_AllBinCounts.txt", sep = "")
s4name <- paste("count_tables/",spec4, "_AllBinCounts.txt", sep = "")
s5name <- paste("count_tables/",spec5, "_AllBinCounts.txt", sep = "")


s1 <- read.table(s1name, header = TRUE)
s2 <- read.table(s2name, header = TRUE)
s3 <- read.table(s3name, header = TRUE)
s4 <- read.table(s4name, header = TRUE)
s5 <- read.table(s5name, header = TRUE)


slist <- list(s1,s2,s3,s4,s5)


for(i in seq(along=slist)){
  count <- slist[[i]]
  count <- count[count$Known >= bin.size,]
  count$GC <- count$GC/count$Known*100
  count[,5:(length(count)-2)] <- ((count[,5:(length(count)-2)]/count$Known) * mb) 
  if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
  if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
  if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
  if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
  if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
  if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}
  
  #count <- count[,!(colnames(count) == "Known")]
  colnames(count)[1:4] <- c("chr", "binID", "start", "end")
  count$binID <- 1:dim(count)[1]
  slist[[i]] <- count
}

KnownS1 <- data.frame(slist[[1]]$binID, slist[[1]]$Known)
KnownS2	<- data.frame(slist[[2]]$binID,	slist[[2]]$Known)
KnownS3 <- data.frame(slist[[3]]$binID, slist[[3]]$Known)
KnownS4  <- data.frame(slist[[4]]$binID,	slist[[4]]$Known)
KnownS5 <- data.frame(slist[[5]]$binID, slist[[5]]$Known)

s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]
s2 <- slist[[2]][,!(colnames(slist[[2]]) == "Known")]
s3 <- slist[[3]][,!(colnames(slist[[3]]) == "Known")]
s4 <- slist[[4]][,!(colnames(slist[[4]]) == "Known")]
s5 <- slist[[5]][,!(colnames(slist[[5]]) == "Known")]

if(rem.un == "yes"){
  if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
  if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
  if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
}
if(rem.un == "yes"){
  if(length(grep("U", s2$chr)) > 0){s2 <- s2[-(grep( "U", s2$chr)),]}
  if(length(grep("_", s2$chr)) > 0){s2 <- s2[-(grep("_", s2$chr)),]}
  if(length(grep("M", s2$chr)) > 0){s2 <- s2[-(grep("M", s2$chr)),]}
}
if(rem.un == "yes"){
  if(length(grep("U", s3$chr)) > 0){s3 <- s3[-(grep( "U", s3$chr)),]}
  if(length(grep("_", s3$chr)) > 0){s3 <- s3[-(grep("_", s3$chr)),]}
  if(length(grep("M", s3$chr)) > 0){s3 <- s3[-(grep("M", s3$chr)),]}
}
if(rem.un == "yes"){
  if(length(grep("U", s4$chr)) > 0){s4 <- s4[-(grep( "U", s4$chr)),]}
  if(length(grep("_", s4$chr)) > 0){s4 <- s4[-(grep("_", s4$chr)),]}
  if(length(grep("M", s4$chr)) > 0){s4 <- s4[-(grep("M", s4$chr)),]}
}
if(rem.un == "yes"){
  if(length(grep("U", s5$chr)) > 0){s5 <- s5[-(grep( "U", s5$chr)),]}
  if(length(grep("_", s5$chr)) > 0){s5 <- s5[-(grep("_", s5$chr)),]}
  if(length(grep("M", s5$chr)) > 0){s5 <- s5[-(grep("M", s5$chr)),]}
}

human <- s1
horse <- s2
bovine <- s3
dog <- s4
mouse <- s5


keep.human <- KnownS1[KnownS1[,2] > keep.limit,1]
keep.horse <- KnownS1[KnownS2[,2] > keep.limit,1]
keep.bovine <- KnownS1[KnownS3[,2] > keep.limit,1]
keep.dog <- KnownS1[KnownS4[,2] > keep.limit,1]
keep.mouse <- KnownS2[KnownS5[,2] > keep.limit,1]

human <- human[human$binID %in% keep.human,]
horse <- horse[horse$binID %in% keep.horse,]
bovine <- bovine[bovine$binID %in% keep.bovine,]
dog <- dog[dog$binID %in% keep.dog,]
mouse <- mouse[mouse$binID %in% keep.mouse,]





pca.human <- prcomp(sqrt(human[,5:ncol(human)]), scale.=T, center=T)


pca.horse <- prcomp(sqrt(horse[,5:ncol(horse)]), scale.=T, center=T)


pca.bovine <- prcomp(sqrt(bovine[,5:ncol(bovine)]), scale.=T, center=T)


pca.dog <- prcomp(sqrt(dog[,5:ncol(dog)]), scale.=T, center=T)


pca.mouse <- prcomp(sqrt(mouse[,5:ncol(mouse)]), scale.=T, center=T)


### fix the above code


#### but can we chnage the colours of our arrows ? so we can highlight the important principal components



# now we can write whatever we want to set the colors. 

# we can make the other names smaller

ycols_human <- rep("black", nrow(pca.human$rotation))
ycols_human[rownames(pca.human$rotation) == "SINE2_MIR"] = "blue"
ycols_human[rownames(pca.human$rotation) == "LINE_L2"] = "blue"
ycols_human[rownames(pca.human$rotation) == "LINE_L1"] = "red"
ycols_human[grep("SINE1_7SL",rownames(pca.human$rotation))] = "darkgreen"
ycols_human[rownames(pca.human$rotation) == "NGenes"] = "purple"

y_cex_human = rep(.9, nrow(pca.human$rotation))
y_cex_human[ycols_human == "black"] =.7



ycols_horse <- rep("black", nrow(pca.horse$rotation))
ycols_horse[rownames(pca.horse$rotation) == "SINE2_MIR"] = "blue"
ycols_horse[rownames(pca.horse$rotation) == "LINE_L2"] = "blue"
ycols_horse[rownames(pca.horse$rotation) == "LINE_L1"] = "red"
ycols_horse[grep("ERE",rownames(pca.horse$rotation))] = "darkgreen"
ycols_horse[rownames(pca.horse$rotation) == "NGenes"] = "purple"

y_cex_horse = rep(.9, nrow(pca.horse$rotation))
y_cex_horse[ycols_horse == "black"] =.7




ycols_bovine <- rep("black", nrow(pca.bovine$rotation))
ycols_bovine[rownames(pca.bovine$rotation) == "SINE2_MIR"] = "blue"
ycols_bovine[rownames(pca.bovine$rotation) == "LINE_L2"] = "blue"
ycols_bovine[rownames(pca.bovine$rotation) == "LINE_L1"] = "red"
ycols_bovine[grep("BOV",rownames(pca.bovine$rotation))] = "darkgreen"
ycols_bovine[grep("ov",rownames(pca.bovine$rotation))] = "darkred"
ycols_bovine[rownames(pca.bovine$rotation) == "NGenes"] = "purple"

y_cex_bovine = rep(.9, nrow(pca.bovine$rotation))
y_cex_bovine[ycols_bovine == "black"] =.7




ycols_dog <- rep("black", nrow(pca.dog$rotation))
ycols_dog[rownames(pca.dog$rotation) == "SINE2_MIR"] = "blue"
ycols_dog[rownames(pca.dog$rotation) == "LINE_L2"] = "blue"
ycols_dog[rownames(pca.dog$rotation) == "LINE_L1"] = "red"
ycols_dog[grep("Can",rownames(pca.dog$rotation))] = "darkgreen"
ycols_dog[rownames(pca.dog$rotation) == "NGenes"] = "purple"

y_cex_dog = rep(.9, nrow(pca.dog$rotation))
y_cex_dog[ycols_dog == "black"] =.7





ycols_mouse <- rep("black", nrow(pca.mouse$rotation))
ycols_mouse[rownames(pca.mouse$rotation) == "SINE2_MIR"] = "blue"
ycols_mouse[rownames(pca.mouse$rotation) == "LINE_L2"] = "blue"
ycols_mouse[rownames(pca.mouse$rotation) == "LINE_L1"] = "red"
ycols_mouse[rownames(pca.mouse$rotation) == "NGenes"] = "purple"
ycols_mouse[grep("7SL",rownames(pca.mouse$rotation))] = "darkgreen"


y_cex_mouse = rep(.9, nrow(pca.mouse$rotation))
y_cex_mouse[ycols_mouse == "black"] =.7




pdf(file = "~/Desktop/GSA_poster/plots_pca.pdf", onefile=T)

reuben.biplot(pca.human$x[,c(1,2)], 
              pca.human$rotation[,c(1,2)], 
              x.col= "grey", 
              y.col=ycols_human, 
              text.col=ycols_human, 
              text.cex=y_cex_human, 
              arrow.lwd=y_cex_human*2.5)
legend("bottomright",legend=spec1, bty="n", cex=1.5)


# green needs to go up
df.x <- pca.horse$x[,c(1,2)]
df.y <- pca.horse$rotation[,c(1,2)]
df.x[,2] = df.x[,2] * -1
df.y[,2] = df.y[,2] * -1
reuben.biplot(df.x[,c(1,2)], 
              df.y[,c(1,2)], 
              x.col= "grey", 
              y.col=ycols_horse, 
              text.col=ycols_horse, 
              text.cex=y_cex_horse, 
              arrow.lwd=y_cex_horse*2.5)
legend("bottomright",legend=spec2, bty="n", cex=1.5)


# blue needs to turn around
df.x <- pca.bovine$x[,c(1,2)] * -1
df.y <- pca.bovine$rotation[,c(1,2)] * -1


reuben.biplot(df.x[,c(1,2)], 
              df.y[,c(1,2)], 
              x.col= "grey", 
              y.col=ycols_bovine, 
              text.col=ycols_bovine, 
              text.cex=y_cex_bovine, 
              arrow.lwd=y_cex_bovine*2.5)
legend("bottomright",legend=spec3, bty="n", cex=1.5)



reuben.biplot(pca.dog$x[,c(1,2)], 
              pca.dog$rotation[,c(1,2)], 
              x.col= "grey", 
              y.col=ycols_dog, 
              text.col=ycols_dog, 
              text.cex=y_cex_dog, 
              arrow.lwd=y_cex_dog*2.5)
legend("bottomright",legend=spec4, bty="n", cex=1.5)



#green needs to go up 
df.x <- pca.mouse$x[,c(1,2)]*-1
df.y <- pca.mouse$rotation[,c(1,2)] *-1

reuben.biplot(df.x[,c(1,2)], 
              df.y[,c(1,2)], 
              x.col= "grey", 
              y.col=ycols_mouse, 
              text.col=ycols_mouse, 
              text.cex=y_cex_mouse, 
              arrow.lwd=y_cex_mouse*2.5)
legend("bottomright",legend=spec5, bty="n", cex=1.5)


dev.off()
# write it that i can get every pca image i need 

# with the correct features pointed out 
pdf(file = "~/Desktop/GSA_poster/legend.pdf")
plot(NA)
legend("bottomleft", c("Ancestral element", "Clade specific SINE", "Number of genes", "L1", "BovB" ), 
       fill = c("blue", "darkgreen", "purple", "red", "darkred"),
       cex=2, box.col = "white")
dev.off()

#####  so lets turn them te right way around and colour them in



# so these are the pca plots of the names of things 
pca.2 <- prcomp(t(scale(dog[,5:ncol(dog)])), scale.=TRUE)
plot(pca.2$x, cex = .5)
text(pca.2$x, labels=rownames(pca.2$x))





my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 19)
my_p_val_palette <- colorRampPalette(c( "white", "black"))(n = 19)
pv <- seq(1,0, by = -.1)
pv.text <- rep("", length(pv))
pv.text[as.integer(quantile(1:(length(pv))))] <- as.character(pv[as.integer(quantile(1:(length(pv))))])
# now i can get the row side colours to give me a p,val


heatmap.2(matrix(seq(0,1, length.out = 100),nrow=10),scale="none", col = my_p_val_palette,dendrogram="n",keysize=2,
          #          symkey = TRUE,
          breaks = seq(from = 0, to = 1, length = 20),
          trace = "none",
          #    symbreaks = TRUE,
          density.info = "none",
          key.title = "",
          key.xlab = "P-value"
)

heatmap.2(matrix(seq(-1.01,1.01, length.out = 100),nrow=10),scale="none", col = my_palette,dendrogram="n",keysize=2,
          #          symkey = TRUE,
          breaks = seq(from = -1, to = 1, length = 20),
          trace = "none",
          #    symbreaks = TRUE,
          density.info = "none",
          key.title = "",
          key.xlab = "Pearson's r"
)



pdf(file = "~/Desktop/GSA_poster/keys.pdf", onefile = T)

layout(matrix(1:3,ncol=3), width = c(1,1,1),height = c(1,1,1))


legend_image <- as.raster(my_palette[19:1])
plot(c(0,2),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Pearson's r", cex.main=2)
text(x=1.5, y = c(1,0.5,0,-0.5, -1), label = c(1,0.5,0,-0.5, -1), cex = 2)
rasterImage(legend_image, 0, -1, 1,1)

plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))


legend_image <- as.raster(my_p_val_palette[19:1])
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "P_value", cex.main=2)
text(x=1.5, y = c(0,0.2,.4,.6,.8, 1), label = c(0,0.2,.4,.6,.8, 1), cex = 2)
rasterImage(legend_image, 0, 0, 1,1)

dev.off()

# I think this script has what I'm looking for 
# hoefully for each of these species we have already converted them over. 


