### sp the next step is to get all the repeat stat figs sorted



rm(list = ls())
setwd("~/Desktop/Domain_manuscript/")
library(GenomicRanges)
library(rtracklayer)


source("Domain_manuscript_scripts/functions.R")
source("Domain_manuscript_scripts/rep_db.R")

spec1 = c("Human", "Chimp", "Rhesus", "Mouse", "Dog")
genome = c("hg19", "panTro4", "rheMac3", "mm9", "canFam3")





pdf(file = "~/Desktop/repDist.pdf", onefile = T, height = 10, width = 12.5)

out <- c(1,3,2,4)
lay = matrix(c(out, out + 4, out +8, out +12, out +16), nrow = 4, byrow = F)
layout(lay)

for(s in 1:length(spec1)){
  rep <- rep_info(spec1 = spec1[s], genome = genome[s])
  repGroup <- repFamStruct(spec1 = spec1[s])
  repGroupName<- unique(repGroup$repType)
  
  for(i in 1:length(repGroupName)){
    teFams <- as.character(repGroup$TEname[repGroup$repType == repGroupName[i]])
    plot(c(0,50), c(0,1), type = "n", main = repGroupName[i], xlab = "percentage mismatch", ylab = "cumulative distribution")
    legend("bottomright", legend = teFams, fill = 1:length(teFams), cex = .7)
    for(te in 1:length(teFams)){
      print(teFams[te])
      te.dens <- rep[[teFams[te]]][,"milliDiv"] / 10
      lines(ecdf(te.dens), col = te)
    }
  }
}
dev.off()



# should we look at element lengths ?
# 

# 



