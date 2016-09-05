### get repets for each species and build dendrograms 


### get comparsions and build heatmaps for that too


### so repeat family stats


### significance test for repeats 
rm(list = ls())


HumanCol <- c(rep("darkgreen",3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))
ChimpCol <- c(rep("darkgreen",3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))
RhesusCol <- c(rep("darkgreen",3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))
MouseCol <- c(rep("darkgreen",4), rep("red", 4), rep("purple", 6), rep("darkblue", 2))
DogCol <- c(rep("darkgreen",5), rep("red", 4), rep("purple", 5), rep("darkblue", 2))


setwd("/Users/labadmin/Desktop/Domain_manuscript")
source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")


load("R_objects/binMaps")


#HumanRef_ChimpQue

head(HumanRef_ChimpQue$binMap)
dim(HumanRef_ChimpQue$HumanRef$data)
dim(HumanRef_DogQue$HumanRef$data)


colFun <- colorRampPalette(c("blue","white" ,"red"))
colFun(20)



pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/selfCompare/Human.pdf", onefile = TRUE)
par(mar = c(5,5,5,5))
reuben.heatmap(cor(HumanRef_ChimpQue$HumanRef$data), scale = "none", col = colFun(20), zlim = c(-1,1), symm = F,
               colsLabs = HumanCol)
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/selfCompare/Chimp.pdf", onefile = TRUE)
par(mar = c(5,5,5,5))
reuben.heatmap(cor(HumanRef_ChimpQue$ChimpQue$data), scale = "none", col = colFun(20), zlim = c(-1,1),  symm = F,
               colsLabs = ChimpCol)
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/selfCompare/Rhesus.pdf", onefile = TRUE)
par(mar = c(5,5,5,5))
reuben.heatmap(cor(HumanRef_RhesusQue$RhesusQue$data), scale = "none", col = colFun(20), zlim = c(-1,1),  symm = F,
               colsLabs = RhesusCol)
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/selfCompare/Mouse.pdf", onefile = TRUE)
par(mar = c(5,5,5,5))
reuben.heatmap(cor(HumanRef_MouseQue$MouseQue$data), scale = "none", col = colFun(20), zlim = c(-1,1),  symm = F,
               colsLabs = MouseCol)
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/selfCompare/Dog.pdf", onefile = TRUE)
par(mar = c(5,5,5,5))
reuben.heatmap(cor(HumanRef_DogQue$DogQue$data), scale = "none", col = colFun(20), zlim = c(-1,1),  symm = F,
               colsLabs = DogCol)
dev.off()

par(mar=c(10,14,10,14))
plot(0,type = "n", axes = FALSE, xlab = "", ylab = "")
image(matrix(seq(-1,1,length.out = 20), nrow=1), col = colFun(20), axes = F)
axis(side = 2, at = seq(0,1,length.out = 5), labels = seq(-1,1,length.out = 5), las = 2)
dev.off()



head(HumanRef_ChimpQue$binMap)
HumanRef_ChimpQue$binMap[1:100,]

# we have the data we just have to do the reordering 

# aggregate 
# and merge 



# we will also need to compare the distributions 
refSpec = "Human"
queSpec = "Chimp"
DataSet <- get(paste(refSpec, "Ref_", queSpec,"Que", sep = ""))

dir.create(paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/", refSpec,queSpec,sep=""))

for(minFrac in c(0,.1,.2,.3,.4,.5)){
  rowChoice <- (1:nrow(DataSet$binMap))[DataSet$binMap$refFrac >= minFrac & DataSet$binMap$queFrac >= minFrac]
  aggQue <- aggregate(DataSet[[paste(queSpec, "Que", sep = "")]]$data[DataSet$binMap$queNo,][rowChoice,] * DataSet$binMap$queFrac[rowChoice], 
                      by = list(DataSet$binMap$refNo[rowChoice]), FUN = sum)
  
  aggHuman <- aggregate(DataSet$binMap$refFrac[rowChoice] ,
                        by = list(DataSet$binMap$refNo[rowChoice]), FUN = sum)
  
  KSpQue <- NULL
  KSdisQue <- NULL
  for(i in 1:(ncol(aggQue)-1)){
    KS <- ks.test(aggQue[,2:ncol(aggQue)][,i]/aggHuman$x,DataSet[[paste(queSpec, "Que", sep = "")]]$data[,i][unique(DataSet$binMap$queNo[rowChoice])])
    KSpQue <- c(KSpQue, KS$p.value)
    KSdisQue <- c(KSdisQue, KS$statistic)
  }
  
  if(queSpec == "Chimp"){
    pdf(file=paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/minMapping/",queSpec,minFrac,".pdf",sep = ""), height = 5,width = 5)
    plot(density((aggQue$L1ME/aggHuman$x)/10000), lwd = 3, xlab = "L1ME density",main="", ylab = "genome segement density")
    lines(density((DataSet[[paste(queSpec, "Que", sep = "")]]$data[,"L1ME"][unique(DataSet$binMap$queNo[rowChoice])])/10000), col = 2, lwd = 3)
    legend("topright",legend = c(paste("min map fraction =",minFrac), paste("P-value =",round(KSpQue[names(aggQue[,2:ncol(aggQue)]) == "L1ME"],digits = 3))),bty = "n")
    dev.off()
  }
  
  KSpRef <- NULL
  KSdisRef <- NULL
  for(i in 1:ncol(DataSet$HumanRef$data)){
    KS <- ks.test(DataSet$HumanRef$data[,i],DataSet$HumanRef$data[,i][unique(DataSet$binMap$refNo[rowChoice])])
    KSpRef <- c(KSpRef, KS$p.value)
    KSdisRef <- c(KSdisRef, KS$statistic)
  }
  
  RefBinNo <- nrow(DataSet$HumanRef$data[unique(DataSet$binMap$refNo[rowChoice]),])/nrow(DataSet$HumanRef$binInfo)
  QueBinNo <- nrow(DataSet[[paste(queSpec, "Que", sep = "")]]$data[unique(DataSet$binMap$queNo[rowChoice]),])/nrow(DataSet[[paste(queSpec, "Que", sep = "")]]$binInfo)
  
  pdf(file = paste("writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/",refSpec,queSpec,"/", refSpec,queSpec,minFrac*10,"frac.pdf", sep = ""))

  #par(mar = c(10,10,10,10))
  heatmap.colours(cor(aggQue[,2:ncol(aggQue)]/aggHuman$x, DataSet$HumanRef$data[aggHuman$Group.1,]), margins = c(10,10),
          scale = "none",col = colFun(20), zlim = c(-1,1),
          RowSideColors = c("white",topo.colors(19))[cut((KSpQue),breaks = c(0,seq(0.05,1,length.out = 20)))],
          ColSideColors = c("white",topo.colors(19))[cut((KSpRef),breaks = c(0,seq(0.05,1,length.out = 20)))],
          xlab = "", ylab = "",RowcolsLabs = get(paste(queSpec,"Col",sep ="")),
          ColcolsLabs = HumanCol
  )
  mtext(side = 1,refSpec,line = 3, cex = 2)
  mtext(side = 4,queSpec,line = 1, cex = 2)
  
  mtext(side = 3,line = -4, text = paste(refSpec,":\n", round(RefBinNo*100,2), "%", sep = "") ,outer = TRUE, at = .00,adj = 0,cex = 1.5)
  mtext(side = 3,line = -4, text = paste(queSpec,":\n", round(QueBinNo*100,2), "%", sep = "") ,outer = TRUE, at = .15,adj = 0,cex = 1.5)
  mtext(side = 3,line = -7, text = paste("Min map fraction:","\n", minFrac*100, "%", sep = "") ,outer = TRUE, at = .00,adj = 0,cex = 1.5)
  dev.off()
}
pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/minMapping/legend.pdf", width = 5,height = 5)
plot.new()
legend("center",legend = c("chimpanzee", "humanised"), fill = c(1,2),title = "Genome", bty = "n", cex = 2)
dev.off()

# species name each level
# proportion of genome copied over at each stage
# this way we can get an accurate profile of TE families and their correlatiojns in other speceis 



pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/legends/legendPval.pdf")
par(mfrow=c(1,2), mar = c(10,8,5,5))
image(matrix(c(seq(0,1,length.out = 1000)), nrow = 1), col = c("white",topo.colors(19)),
      xaxt = "n", yaxt = "n",ylim = c(0,1), main = "P-value", cex.main = 2.5)
box(which = "plot", lwd = 5)
axis(2,seq(0,1,length.out = 21)[c(1,6,11,16,21)], labels = FALSE, las = 2, lwd.ticks = 5)
mtext(text = round(seq(0,1,length.out = 21),2)[c(1,6,11,16,21)],side = 2,line = 2,
      at = seq(0,1,length.out = 21)[c(1,6,11,16,21)], las = 2, cex = 2)
axis(4,seq(0,1,length.out = 21)[2], labels = FALSE, las = 2, col = 2,lwd = 5)
mtext(text = round(seq(0,1,length.out = 21),2)[2],side = 4,line = 2,
      at = seq(0,1,length.out = 21)[2], las = 2, cex = 2)
abline(h=.05, col = 2, lwd = 5)
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/legends/legendCor.pdf")
par(mfrow=c(1,2), mar = c(10,8,5,5))
image(matrix(seq(-1,1,length.out = 20), nrow = 1), col = colFun(20), ylim = c(0,1), 
      xaxt = "n", yaxt = "n", main = "cor", cex.main = 2.5)
box(lwd = 5)
axis(2,seq(0,1,.25),labels = FALSE,las = 2, lwd= 5)
mtext(side = 2,at = seq(0,1,.25),text  = seq(-1,1,length.out = 5),las = 2,line = 2,cex = 2)
dev.off()


pdf(file = "writing/round2_20160503/draftsTex/supmaterial/TexFigs/supFig/legends/groups.pdf")
layout(1)
plot(1,type = "n", axes = FALSE,xlab = "", ylab = "")
legend("center", legend = c("new SINE", "new L1", "old L1", "ancient"), title = "retrotransposon\ngroups",
       fill = c("darkgreen", "purple", "red", "darkblue"), box.lwd = -1, cex = 2,box.col = "white")
dev.off()



# so we got 
# 1. raw distribution of que
# 2. raw distribution of que subset
# 3. remodeled distribution of que (is a subset)
# 4. raw distribution of ref
# 5. raw distribution of ref (subset)





### comparisons
## 1, 2
## 1, 3
## 2, 3
## 4, 5





# so how best to get this data
# we could go just to grey scale at heatmaps

### get the colours correct 


par(mar=c(5,5,5,5))

pcaChimp <- prcomp(aggQue[,2:ncol(aggQue)]/aggHuman$x, scale. = T)
reuben.biplot(x = pcaChimp$x[,1:2],y = pcaChimp$rotation[,1:2] ,text.col = DogCol,y.col = DogCol)


## maybe we could colour points based on their mapping

# get the genome conservation stuff working 
# then we will be close to done for our RTN comparisons. 

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/human/chimp.pdf",width = 5,height = 5)
# Chimp
hFragAggChimp <- aggregate(HumanRef_ChimpQue$binMap$refFrac, by = list(HumanRef_ChimpQue$binMap$refNo), sum)
cuts <- cut(hFragAggChimp$x,breaks = seq(0,.7,length.out = 20))
pca <- prcomp(HumanRef_ChimpQue$HumanRef$data[hFragAggChimp$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/human/rhesus.pdf",width = 5,height = 5)
# Rhesus
hFragAggRhesus <- aggregate(HumanRef_RhesusQue$binMap$refFrac, by = list(HumanRef_RhesusQue$binMap$refNo), sum)
cuts <- cut(hFragAggRhesus$x,breaks = seq(0,.7,length.out = 20))
pca <- prcomp(HumanRef_RhesusQue$HumanRef$data[hFragAggRhesus$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/human/mouse.pdf",width = 5,height = 5)
# Mouse
hFragAggMouse <- aggregate(HumanRef_MouseQue$binMap$refFrac, by = list(HumanRef_MouseQue$binMap$refNo), sum)
cuts <- cut(hFragAggMouse$x,breaks = seq(0,.7,length.out = 20))
pca <- prcomp(HumanRef_MouseQue$HumanRef$data[hFragAggMouse$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/human/dog.pdf",width = 5,height = 5)
# Dog
hFragAggDog <- aggregate(HumanRef_DogQue$binMap$refFrac, by = list(HumanRef_DogQue$binMap$refNo), sum)
cuts <- cut(hFragAggDog$x,breaks = 20)
pca <- prcomp(HumanRef_DogQue$HumanRef$data[hFragAggDog$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()






pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/legends/legendMapping.pdf")
par(mfrow=c(1,2), mar = c(10,8,5,5))
image(matrix(seq(-1,1,length.out = 20), nrow = 1), col = colFun(20), ylim = c(0,1), 
      xaxt = "n", yaxt = "n", main = "mapping\nfraction", cex.main = 2.5)
box(lwd = 5)
axis(2,seq(0,1,.25),labels = FALSE,las = 2, lwd= 5)
mtext(side = 2,at = seq(0,1,.5),text  = seq(0,.7,length.out = 3),las = 2,line = 2,cex = 2)
dev.off()

# start to build up the sup a bit more
# now that we have the species conservation scores we need. 

### we need to look at the other direction
layout(1)


pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/nonHuman/chimp.pdf",width = 5,height = 5)
cFragAggHuman <- aggregate(HumanRef_ChimpQue$binMap$queFrac, by = list(HumanRef_ChimpQue$binMap$queNo), sum)
cuts <- cut(10^(cFragAggHuman$x),breaks = 20)
pca <- prcomp(HumanRef_ChimpQue$ChimpQue$data[cFragAggHuman$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
#reuben.biplot(data.frame(pca$x[,2]*1,pca$x[,1]*1), data.frame(pca$rotation[,2]*1,pca$rotation[,1]*1), x.col=colFun(20)[cuts])
plot(data.frame(pca$x[,2]*1,pca$x[,1]*1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/nonHuman/rhesus.pdf",width = 5,height = 5)
rFragAggHuman <- aggregate(HumanRef_RhesusQue$binMap$queFrac, by = list(HumanRef_RhesusQue$binMap$queNo), sum)
cuts <- cut(10^(rFragAggHuman$x),breaks = 20)
pca <- prcomp(HumanRef_RhesusQue$RhesusQue$data[rFragAggHuman$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
#reuben.biplot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), data.frame(pca$rotation[,2]*-1,pca$rotation[,1]*-1), x.col=colFun(20)[cuts])
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/nonHuman/mouse.pdf",width = 5,height = 5)
mFragAggHuman <- aggregate(HumanRef_MouseQue$binMap$queFrac, by = list(HumanRef_MouseQue$binMap$queNo), sum)
cuts <- cut(mFragAggHuman$x,breaks = seq(0,.7,length.out = 20))
pca <- prcomp(HumanRef_MouseQue$MouseQue$data[mFragAggHuman$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
#reuben.biplot(data.frame(pca$x[,2]*-1,pca$x[,1]*1), data.frame(pca$rotation[,2]*-1,pca$rotation[,1]*1), x.col=colFun(20)[cuts])
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/PCA/nonHuman/dog.pdf",width = 5,height = 5)
dFragAggHuman <- aggregate(HumanRef_DogQue$binMap$queFrac, by = list(HumanRef_DogQue$binMap$queNo), sum)
cuts <- cut(dFragAggHuman$x,breaks = seq(0,.7,length.out = 20))
pca <- prcomp(HumanRef_DogQue$DogQue$data[dFragAggHuman$Group.1,], scale. = T)
par(mar=c(5,5,5,5))
#reuben.biplot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), data.frame(pca$rotation[,2]*-1,pca$rotation[,1]*-1), x.col=colFun(20)[cuts])
plot(data.frame(pca$x[,2]*-1,pca$x[,1]*-1), col = colFun(20)[cuts], pch = 16, cex = .5,
     xlab = "Ancient PC", ylab ="new SINE PC")
dev.off()



###### barplots 

43076201+
  28322032+
  18534030+
  15557796


sum(34095135,
    27706210,
    13378123,
    10387383)


# repeat families from each species 
HumanPercent <- (colSums(HumanRef_ChimpQue$HumanRef$data)*100) /sum(HumanRef_ChimpQue$HumanRef$binInfo$Known)
ChimpPercent <- (colSums(HumanRef_ChimpQue$ChimpQue$data)*100) /sum(HumanRef_ChimpQue$ChimpQue$binInfo$Known)
RhesusPercent <- (colSums(HumanRef_RhesusQue$RhesusQue$data)*100) /sum(HumanRef_RhesusQue$RhesusQue$binInfo$Known)
MousePercent <- (colSums(HumanRef_MouseQue$MouseQue$data)*100) /sum(HumanRef_MouseQue$MouseQue$binInfo$Known)
DogPercent <- (colSums(HumanRef_DogQue$DogQue$data)*100) /sum(HumanRef_DogQue$DogQue$binInfo$Known)

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/Human.pdf", width = 10,height = 3)
barplot(rep(NA,length(HumanPercent)), col = HumanCol,las = 2, ylim= c(0,7),axes = F,names.arg = "")
grid()
barplot(HumanPercent*100, col = HumanCol,las = 2, ylim= c(0,7), ylab = c("genome coverage (%)"),axes = F,plot = T, add = T)
axis(side = 2,at = (0:7), labels = c(0,NA,2,NA,4,NA,6,NA))
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/Chimp.pdf", width = 10,height = 3)
barplot(rep(NA,length(ChimpPercent)), col = HumanCol,las = 2, ylim= c(0,7),axes = F,names.arg = "")
grid()
barplot(ChimpPercent* 100, col = ChimpCol,las = 2, ylim= c(0,7), ylab = c("genome coverage (%)"), add = T, axes = F)
axis(side = 2,at = (0:7), labels = c(0,NA,2,NA,4,NA,6,NA))
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/Rhesus.pdf", width = 10,height = 3)
barplot(rep(NA,length(RhesusPercent)), col = HumanCol,las = 2, ylim= c(0,7),axes = F,names.arg = "")
grid()
barplot(RhesusPercent* 100, col = RhesusCol,las = 2, ylim= c(0,7), ylab = c("genome coverage (%)"), add = T, axes = F)
axis(side = 2,at = (0:7), labels = c(0,NA,2,NA,4,NA,6,NA))
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/Mouse.pdf", width = 10,height = 3)
barplot(rep(NA,length(MousePercent)), col = HumanCol,las = 2, ylim= c(0,7),axes = F,names.arg = "")
grid()
barplot(MousePercent* 100, col = MouseCol,las = 2, ylim= c(0,7), ylab = c("genome coverage (%)"), add = T, axes = F)
axis(side = 2,at = (0:7), labels = c(0,NA,2,NA,4,NA,6,NA))
dev.off()

pdf(file = "writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/Dog.pdf", width = 10,height = 3)
barplot(rep(NA,length(DogPercent)), col = HumanCol,las = 2, ylim= c(0,7),axes = F,names.arg = "")
grid()
barplot(DogPercent * 100, col = DogCol,las = 2, ylim= c(0,7), ylab = c("genome coverage (%)"), add = T, axes = F)
axis(side = 2,at = (0:7), labels = c(0,NA,2,NA,4,NA,6,NA))
dev.off()

# repeat groups from each speceis 

newSINE <- c(Human = sum(HumanPercent[HumanCol == "darkgreen"]), 
             Chimp = sum(ChimpPercent[ChimpCol == "darkgreen"]), 
             Rhesus = sum(RhesusPercent[RhesusCol == "darkgreen"]), 
             Mouse = sum(MousePercent[MouseCol == "darkgreen"]), 
             Dog = sum(DogPercent[DogCol == "darkgreen"]))

newLINE <- c(Human = sum(HumanPercent[HumanCol == "purple"]), 
             Chimp = sum(ChimpPercent[ChimpCol == "purple"]), 
             Rhesus = sum(RhesusPercent[RhesusCol == "purple"]), 
             Mouse = sum(MousePercent[MouseCol == "purple"]), 
             Dog = sum(DogPercent[DogCol == "purple"]))

oldLINE <- c(Human = sum(HumanPercent[HumanCol == "red"]), 
             Chimp = sum(ChimpPercent[ChimpCol == "red"]), 
             Rhesus = sum(RhesusPercent[RhesusCol == "red"]), 
             Mouse = sum(MousePercent[MouseCol == "red"]), 
             Dog = sum(DogPercent[DogCol == "red"]))

ancient <- c(Human = sum(HumanPercent[HumanCol == "darkblue"]), 
             Chimp = sum(ChimpPercent[ChimpCol == "darkblue"]), 
             Rhesus = sum(RhesusPercent[RhesusCol == "darkblue"]), 
             Mouse = sum(MousePercent[MouseCol == "darkblue"]), 
             Dog = sum(DogPercent[DogCol == "darkblue"]))


pdf(file="writing/round2_20160503/draftsTex/TexFigs/supFig/GenomeContent/all.pdf")
barplot(rep(NA,20), las = 2, space = c(c(rep(0.1,5),.4),  c(rep(0.1,4),.4), c(rep(0.1,4),.4), c(rep(0.1,4))), axes = F,
        col = c(rep("darkgreen",5), rep("purple",5), rep("red",5), rep("darkblue",5)), ylab = "", ylim = c(0,20))
grid()
barplot(c(newSINE,newLINE,oldLINE,ancient) *100, las = 2, space = c(c(rep(0.1,5),.4),  c(rep(0.1,4),.4), c(rep(0.1,4),.4), c(rep(0.1,4))),
        col = c(rep("darkgreen",5), rep("purple",5), rep("red",5), rep("darkblue",5)), ylab = "genome coverage (%)", ylim = c(0,20), add = T)
dev.off()








