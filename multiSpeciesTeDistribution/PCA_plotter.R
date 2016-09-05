# PCA plotter


rm(list = ls())
load("~/Desktop/Domain_manuscript/R_objects/PCA_species")


source(file="~/Desktop/Domain_manuscript/Domain_manuscript_scripts/functions.R")
library(GenomicRanges)
library(rtracklayer)
library(circlize)



pdf(file="~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/PCAspec.pdf", onefile = T, height = 5,width = 5)

layout(matrix(c(1,2,3,4), nrow = 2))
par(mar=c(2,2,2,2))

ycol <- c(rep("darkgreen", 3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))
reuben.biplot(x=ChimpPCA$x,y=ChimpPCA$rotation, cex=.2, arrow.lwd=2,text.cex = .001,
              y.col=ycol, 
              text.col=ycol,
              xlim = c(-7,7),
              ylim = c(-7,7),
              ratio = .06, ylab = "", las = 1,xaxt = "n", yaxt = "n",
#               xlab = paste(as.character(round(ChimpPCA$importance$ancient_PC, digits=2)*100), "%"),
#               ylab = paste(as.character(round(ChimpPCA$importance$new_SINE_PC, digits=2)*100), "%"),
x.col = "grey60"
)
legend("topright", legend="Chimp", cex = 1, bty = "n")
legend("bottomright", legend = c("new SINE", "new L1", "old L1", "ancient"), 
       fill = c("darkgreen", "purple", "red", "darkblue"), bty = "n", cex = .8)
par(new = TRUE)
plot(1, type = "n", xlim = c(-7,7), ylim = c(-7,7), axes = FALSE)
axis(side = 1,c(-5,0,5))
axis(side = 2,c(-5,0,5))





ycol <- c(rep("darkgreen", 3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))
reuben.biplot(x=RhesusPCA$x,y=RhesusPCA$rotation, cex=.2, arrow.lwd=2,text.cex = .001,
              y.col=ycol, 
              text.col=ycol,
              xlim = c(-7,7),
              ylim = c(-7,7),
              ratio = .06, xaxt = "n", yaxt = "n", ylab = "",
#               xlab = paste(as.character(round(RhesusPCA$importance$ancient_PC, digits=2)*100), "%"),
#               ylab = paste(as.character(round(RhesusPCA$importance$new_SINE_PC, digits=2)*100), "%")    
x.col = "grey60"
)
legend("topright", legend="Rhesus", cex = 1, bty = "n")
par(new = TRUE)
plot(1, type = "n", xlim = c(-7,7), ylim = c(-7,7), axes = FALSE)
axis(side = 1,c(-5,0,5))
axis(side = 2,c(-5,0,5))



ycol <- c(rep("darkgreen", 4), rep("red", 4), rep("purple", 6), rep("darkblue", 2))
reuben.biplot(x=cbind(MousePCA$x[,1],MousePCA$x[,2]) ,y=cbind(MousePCA$rotation[,1],MousePCA$rotation[,2]), cex=.2, arrow.lwd=2,text.cex = .001,
              y.col=ycol, 
              text.col=ycol,
              xlim = c(-7,7),
              ylim = c(-7,7),
              ratio = .06, xaxt = "n", yaxt = "n", ylab = "", 
#               xlab = paste(as.character(round(MousePCA$importance$ancient_PC, digits=2)*100), "%"),
#               ylab = paste(as.character(round(MousePCA$importance$new_SINE_PC, digits=2)*100), "%")  
x.col = "grey60"
)
legend("topright", legend="Mouse", cex = 1, bty = "n")
par(new = TRUE)
plot(1, type = "n", xlim = c(-7,7), ylim = c(-7,7), axes = FALSE)
axis(side = 1,c(-5,0,5))
axis(side = 2,c(-5,0,5))


ycol <- c(rep("darkgreen", 5), rep("red", 4), rep("purple", 5), rep("darkblue", 2))
reuben.biplot(x=DogPCA$x,y=DogPCA$rotation, cex=.2, arrow.lwd=2,text.cex = .001,
              y.col=ycol, 
              text.col=ycol,
              xlim = c(-7,7),
              ylim = c(-7,7),
              ratio = .06, xaxt = "n", yaxt = "n", ylab = "",
    #          xlab = paste(as.character(round(DogPCA$importance$ancient_PC, digits=2)*100), "%"),
     #         ylab = paste(as.character(round(DogPCA$importance$new_SINE_PC, digits=2)*100), "%")     
    x.col = "grey60"
              )
legend("topright", legend="Dog", cex = 1, bty = "n")
par(new = TRUE)
plot(1, type = "n", xlim = c(-7,7), ylim = c(-7,7), axes = FALSE)
axis(side = 1,c(-5,0,5))
axis(side = 2,c(-5,0,5))


dev.off()


#pdf(file = "~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/PCAhuman.pdf", height = 5, width = 5)





A <- HumanPCA$binInfo
nos <- (1:nrow(A))[A$chr != "chrY"]


which <- GRanges(seqnames=Rle(HumanPCA$binInfo$chr[nos]), 
                   ranges=IRanges(start = HumanPCA$binInfo$start[nos], end = HumanPCA$binInfo$end[nos])
)
repliTime <- import("~/Desktop/Domain_manuscript/Data/UW_repliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format = "bw", which = which)

ol <- as.matrix(findOverlaps(which, repliTime))
ol.agg <- aggregate(x = elementMetadata(repliTime)$score[ol[,2]], by = list(ol[,1]), FUN=mean)


f = colorRamp2(breaks = c((min((ol.agg$x))), mean((ol.agg$x)), max(abs(ol.agg$x))), colors = c("blue4", "white", "green4"))
ycol <- c(rep("darkgreen", 3), rep("red", 4), rep("purple", 4), rep("darkblue", 2))



pdf(file = "~/Desktop/Domain_manuscript/plots/FIg1/Human_repli_pca.pdf", onefile = T, width = 5, height = 5)

layout(matrix(c(1,2,3,4), nrow = 2))
par(mar=c(2,2,2,2))


reuben.biplot(x=HumanPCA$x[nos,],y=HumanPCA$rotation[nos,], cex=.2, arrow.lwd=2, text.cex = 0.001,
              y.col=ycol, 
              text.col=ycol,
#               xlab = paste(as.character(round(HumanPCA$importance$ancient_PC, digits=2)*100), "%"),
#               ylab = paste(as.character(round(HumanPCA$importance$new_SINE_PC, digits=2)*100), "%"),
      x.col = "grey", xlim = c(-7,7), ylim = c(-7,7),ratio = .06, xaxt = "n",yaxt = "n"
)
legend("topright", legend="Human", cex = 1, bty = "n")
par(new=TRUE)
plot(1,type = "n", axes =FALSE, xlim = c(-7,7), ylim = c(-7,7))
axis(1,c(-5,0,5))
axis(2,c(-5,0,5))



x=HumanPCA$x[nos,]
y=HumanPCA$rotation[nos,]
unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
                                abs(max(x, na.rm = TRUE)))
rangx1 <- unsigned.range(x[, 1L])
rangx2 <- unsigned.range(x[, 2L])
rangy1 <- unsigned.range(y[, 1L])
rangy2 <- unsigned.range(y[, 2L])
  xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)

plot.new()
  
plot(x=HumanPCA$x[nos,1],y=HumanPCA$x[nos,2], cex=.4, pch = 16, xlim=c(-7,7),ylim=c(-7,7),
              xlab = paste(as.character(round(HumanPCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(HumanPCA$importance$new_SINE_PC, digits=2)*100), "%"),
              col = f(ol.agg$x), xaxt = "n", yaxt="n"
)
legend("topright", legend="Replication timing", cex = 1, bty = "n")
axis(1,c(-5,0,5))
axis(2,c(-5,0,5))

dev.off()




pdf(file = "~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/legends.pdf")
layout(c(1,2))
par(mar=c(13,15,3,15))
image(t(matrix(1:5, nrow=1)),col=f( seq(min(ol.agg$x),  max(ol.agg$x)) ),   xaxt = "n", yaxt = "n")
axis(side = 1, at = c(0, 1), labels = c("late", "early"), las=1)



par(mar= c(2,1,1,1))
plot.new()
legend("center", c("old LINE","new LINE", "ancestral", "new SINE"), fill = c("red","purple", "darkblue", "darkgreen"),horiz = F)
legend("bottom", c("old LINE","new LINE", "ancestral", "new SINE"), fill = c("red","purple", "darkblue", "darkgreen"),horiz = T)

dev.off()



# I think mouse needs to be turned around 





