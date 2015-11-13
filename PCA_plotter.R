# PCA plotter


rm(list = ls())
load("~/Desktop/Domain_manuscript/R_objects/PCA_species")


source(file="~/Desktop/element_curves/element_curves_scripts/functions.R")


pdf(file="~/Desktop/Domain_manuscript/plots/FIg1/PCA_wiggle/PCA.pdf")

layout(matrix(c(1,2,3,4), nrow = 2))
par(mar=c(4,4,4,4))

ycol <- c(rep("darkgreen", 5), rep("red", 5), rep("purple", 4), rep("darkblue", 2))
reuben.biplot(x=DogPCA$x,y=DogPCA$rotation, cex=.2, arrow.lwd=2,
              y.col=ycol, 
              text.col=ycol,
              xlab = paste(as.character(round(DogPCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(DogPCA$importance$new_SINE_PC, digits=2)*100), "%")                  
              )
legend("bottomright", legend="Dog", cex = 2, bty = "n")




ycol <- c(rep("darkgreen", 4), rep("red", 5), rep("purple", 5), rep("darkblue", 2))
reuben.biplot(x=MousePCA$x,y=MousePCA$rotation, cex=.2, arrow.lwd=2,
              y.col=ycol, 
              text.col=ycol,
              xlab = paste(as.character(round(MousePCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(MousePCA$importance$new_SINE_PC, digits=2)*100), "%")                  
)
legend("bottomright", legend="Mouse", cex = 2, bty = "n")



ycol <- c(rep("darkgreen", 3), rep("red", 5), rep("purple", 3), rep("darkblue", 2))
reuben.biplot(x=ChimpPCA$x,y=ChimpPCA$rotation, cex=.2, arrow.lwd=2,
              y.col=ycol, 
              text.col=ycol,
              xlab = paste(as.character(round(ChimpPCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(ChimpPCA$importance$new_SINE_PC, digits=2)*100), "%")                  
)
legend("bottomright", legend="Chimp", cex = 2, bty = "n")



ycol <- c(rep("darkgreen", 3), rep("red", 5), rep("purple", 3), rep("darkblue", 2))
reuben.biplot(x=RhesusPCA$x,y=RhesusPCA$rotation, cex=.2, arrow.lwd=2,
              y.col=ycol, 
              text.col=ycol,
              xlab = paste(as.character(round(RhesusPCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(RhesusPCA$importance$new_SINE_PC, digits=2)*100), "%")                  
)
legend("bottomright", legend="Rhesus", cex = 2, bty = "n")


dev.off()



layout(matrix(1))

ycol <- c(rep("darkgreen", 3), rep("red", 5), rep("purple", 3), rep("darkblue", 2))
reuben.biplot(x=HumanPCA$x,y=HumanPCA$rotation, cex=.2, arrow.lwd=2,
              y.col=ycol, 
              text.col=ycol,
              xlab = paste(as.character(round(HumanPCA$importance$ancient_PC, digits=2)*100), "%"),
              ylab = paste(as.character(round(HumanPCA$importance$new_SINE_PC, digits=2)*100), "%")                  
)
legend("bottomright", legend="Human", cex = 2, bty = "n")








