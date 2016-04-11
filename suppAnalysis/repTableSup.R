
rm(list = ls())
genome = "rheMac3"
spec1 = "Rhesus"

setwd("~/Desktop/Domain_manuscript/")


rep_name <- paste( "~/Desktop/topological_domains/TD_Data/rmsk/", genome,"_rmsk.txt" , sep = "")


rep <- read.table(file= rep_name,
                  colClasses= c("integer", "integer", "integer", "integer", "integer",  # some stats
                                "character", "integer", "integer",  # repeat coordinates
                                "integer",  # genome left
                                "character",   # strand
                                "factor", "factor", "factor",     # classifications
                                "integer", "integer", "integer", "integer"),        # more coordinates
                  comment.char="",
                  header = T, fill = T
)


# maybe we have a tree/list object and work our way down to the leaves to get all the stats
# then we can convert it all to a csv file later, since we know how to get the leaves for 
# each species it shouldnt be too hard





MIR <- rep[grep("MIR", rep$repFamily),]
L2 <- rep[grep("L2", rep$repFamily),]
L1 <- rep[grep("L1", rep$repFamily),]
L1ME <- L1[grep("ME", L1$repName),]
L1MD <- L1[grep("MD", L1$repName),]
L1MC <- L1[grep("MC", L1$repName),]
L1MB <- L1[grep("MB", L1$repName),]
L1MA <- L1[grep("MA", L1$repName),]

if(spec1 == "Human"){
  AluJ <- rep[grep("AluJ", rep$repName),]
  AluS <- rep[grep("AluS", rep$repName),]
  AluY <- rep[grep("AluY", rep$repName),]
  L1PB <- L1[grep("PB", L1$repName),]
  L1PA <- L1[grep("PA", L1$repName),]
  L1HS <- L1[grep("HS", L1$repName),] 
}

if(spec1 == "Mouse"){
  ALU <- rep[grep("Alu", rep$repFamily),]
  PB <- ALU[grep("PB", ALU$repName),]
  B1 <- rbind(ALU[grep("B1_", ALU$repName),], ALU[grep("B1F", ALU$repName),])
  B2all <- rep[grep(pattern="B2", rep$repFamily),]
  B2 <- B2all[grep("B2", B2all$repName),]
  B3 <- B2all[grep("B3", B2all$repName),]
  B4all <- rep[grep(pattern="B4", rep$repFamily),]
  B4 <- B4all[-(grep(pattern="RSINE", B4all$repName)),]
  L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),]
  Lx <- L1new[grep("Lx", L1new$repName),]
  L1Md <- L1new[grep("L1Md", L1new$repName),]
  L1_Mus <- L1new[grep("L1_Mus", L1new$repName),]
  L1_Mur <- L1new[grep("L1_Mur", L1new$repName),]
  L1_Mm <- L1new[grep("L1_Mm", L1new$repName),]
}

if(spec1 == "Chimp"){
  Alu <- rep[grep("Alu", rep$repFamily),]
  AluJ <- Alu[grep("AluJ", Alu$repName),]
  AluS <- Alu[grep("AluS", Alu$repName),]
  AluY <- Alu[grep("AluY", Alu$repName),]
  L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),] 
  L1PB <- L1new[grep("PB", L1new$repName),]
  L1PA <- L1new[grep("PA", L1new$repName),]
  L1Pt <- L1new[grep("Pt", L1new$repName),] 
}

if(spec1 == "Rhesus"){
  Alu <- rep[grep("Alu", rep$repFamily),]
  AluJ <- Alu[grep("AluJ", Alu$repName),]
  AluS <- Alu[grep("AluS", Alu$repName),]
  AluY <- Alu[grep("AluY", Alu$repName),]
  L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),] 
  L1PB <- L1new[grep("PB", L1new$repName),]
  L1PA <- L1new[grep("PA", L1new$repName),]
  L1RS <- rbind(L1new[grep("_RS", L1new$repName),], L1new[grep("PREC", L1new$repName),])
}

if(spec1 == "Dog"){
  Lys <- rep[grep("tRNA-Lys", rep$repFamily),]
  SINEC_c <- Lys[grep("SINEC_c", Lys$repName),]
  SINEC_b <- Lys[grep("SINEC_b", Lys$repName),]
  SINEC_a <- Lys[grep("SINEC_a", Lys$repName),]
  SINEC_old <- Lys[grep("SINEC_old", Lys$repName),]
  SINEC_Cf <- Lys[grep("SINEC_Cf", Lys$repName),]
  L1new <- L1[!(as.character(L1$repName) %in% c(as.character(L1ME$repName), as.character(L1MD$repName), as.character(L1MC$repName), as.character(L1MB$repName), as.character(L1MA$repName))  ),]
  L1_Carn <- L1new[grep("L1_Carn", L1new$repName), ]
  L1_Canid <- L1new[grep("L1_Canid", L1new$repName), ]
  L1_Canis <- L1new[grep("L1_Canis", L1new$repName), ]
  L1_Cf <- L1new[grep("L1_Cf", L1new$repName), ]
}



specList<- list(Ancient = list(MIR = NULL, 
                               L2 = NULL),
                old_L1 = list(L1ME = NULL,
                              L1MD = NULL,
                              L1MC = NULL,
                              L1MB = NULL)
                )

if(spec1 == "Human"){
  specList$new_SINE <- list(AluJ = NULL, 
                            AluS = NULL,
                            AluY = NULL)
  specList$new_L1 <- list(L1MA = NULL,
                          L1PB = NULL,
                          L1PA = NULL, 
                          L1HS = NULL)
  
}

if(spec1 == "Chimp"){
  specList$new_SINE <- list(AluJ = NULL, 
                            AluS = NULL,
                            AluY = NULL)
  specList$new_L1 <- list(L1MA = NULL,
                          L1PB = NULL,
                          L1PA = NULL, 
                          L1Pt = NULL)
}

if(spec1 == "Rhesus"){
  specList$new_SINE <- list(AluJ = NULL, 
                            AluS = NULL,
                            AluY = NULL)
  specList$new_L1 <- list(L1MA = NULL,
                          L1PB = NULL,
                          L1PA = NULL, 
                          L1RS = NULL)
}

if(spec1 == "Dog"){
  specList$new_SINE <- list(SINEC_c = NULL, 
                            SINEC_b = NULL,
                            SINEC_a = NULL,
                            SINEC_old = NULL,
                            SINEC_Cf = NULL)
  specList$new_L1 <- list(L1MA = NULL,
                          L1_Carn = NULL,
                          L1_Canid = NULL, 
                          L1_Canis = NULL,
                          L1_Cf = NULL)
}

if(spec1 == "Mouse"){
  specList$new_SINE <- list(PB = NULL, 
                            B1 = NULL,
                            B2 = NULL,
                            B3 = NULL,
                            B4 = NULL)
  specList$new_L1 <- list(L1MA = NULL,
                          Lx = NULL,
                          L1Md = NULL, 
                          L1_Mus = NULL,
                          L1_Mur = NULL,
                          L1_Mm = NULL)
}


### as we run through here lets also start to build a dataframe we can convert to a csv later on for human
specDF <- NULL

for(l1 in 1:length(names(specList))){
  teL2 <- names(specList[[names(specList)[l1]]])
  for(l3 in 1:length(teL2)){
    Te <- get(teL2[l3])
    teSum <- summary(Te$repName)[summary(Te$repName)>0]
    for(l4 in 1:length(teSum)){
      coverage <- Te[Te$repName == names(teSum)[l4], c("genoStart", "genoEnd")]
      coverage= sum(coverage$genoEnd - coverage$genoStart + 1)
      freq = teSum[l4]
      names(freq) = NULL
      specList[[l1]][[teL2[l3]]][[names(teSum[l4])]] <- data.frame(frequency = freq, coverage = coverage)
      
      specDF <- rbind(specDF, c(repGroup = names(specList)[l1],repFam =  teL2[l3], repID = names(teSum[l4]), frequency = freq, coverage = coverage))
      
    }
  }
}
specDF <- as.data.frame(specDF)



specDF$repGroup = as.character(specDF$repGroup)
specDF[duplicated(specDF$repGroup),"repGroup"] = ""

specDF$repFam = as.character(specDF$repFam)
specDF[duplicated(specDF$repFam),"repFam"] = ""

write.csv(x = specDF, file = paste("supTable/", spec1,"_", genome,"_repTab.csv", sep = ""))
          
          
          
