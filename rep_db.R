# species info 

# write a function which can decleare the variables needed for the analysis 
# tell it which species to get the info, it will make a bunch of objects 
# these will be available via TE.names


rep_info <- function(spec1, genome){
  
  rep_name <- paste( "~/Desktop/Domain_manuscript/Data/rmskTable/", genome,"_rmsk.txt" , sep = "")
  
  
  rep <- read.table(file= rep_name,
                    colClasses= c("integer", "integer", "integer", "integer", "integer",  # some stats
                                  "character", "integer", "integer",  # repeat coordinates
                                  "integer",  # genome left
                                  "character",   # strand
                                  "factor", "factor", "factor",     # classifications
                                  "integer", "integer", "integer", "integer"),        # more coordinates
                    col.names= c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                                 "genoName", "genoStart", "genoEnd", "genoLeft",
                                 "strand", "repName", "repClass", "repFamily", 
                                 "repStart", "repEnd", "repLeft", "id"),
                    comment.char="#",
                    header = FALSE, fill = TRUE
  )
  
  
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
  
  
  # just get some more TE clases over the genome so we can get a more refined look at whats happening
  if(spec1 == "Human"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1HS", "L2", "MIR")
  }
  if(spec1 == "Mouse"){
    TE.names <- c("B1", "B2", "B3", "B4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "Lx", "L1_Mur", "L1_Mus","L1Md" ,"L1_Mm", "L2", "MIR")
  }
  if(spec1 == "Chimp"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1Pt", "L2", "MIR")
  }
  if(spec1 == "Rhesus"){
    TE.names <- c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1RS", "L2", "MIR")
  }
  if(spec1 == "Dog"){
    TE.names <- c("SINEC_old", "SINEC_c", "SINEC_b", "SINEC_a", "SINEC_Cf", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_Carn", "L1_Canid", "L1_Canis", "L1_Cf", "L2", "MIR")
  }
  
  
  sorted.TE <- NULL
  for(i in TE.names){
    sorted.TE = c(sorted.TE, list(get(i)))
  }
  names(sorted.TE) <- TE.names
  
  return(sorted.TE)
  
}

#### a function where we can get tables that describe the family structure for each TE

repFamStruct <- function(spec1){
  if(spec1 == "Human"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1HS", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Mouse"){
    repFamStruct <- data.frame(TEname = c("B1", "B2", "B3", "B4", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "Lx", "L1_Mur", "L1_Mus","L1Md" ,"L1_Mm", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Chimp"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1Pt", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
    )
  }
  if(spec1 == "Rhesus"){
    repFamStruct <- data.frame(TEname = c("AluJ", "AluS", "AluY", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1PB", "L1PA", "L1RS", "L2", "MIR"),
                               repType = c("new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")    
                               )
  }
  if(spec1 == "Dog"){
    repFamStruct <- data.frame(TEname = c("SINEC_old", "SINEC_c", "SINEC_b", "SINEC_a", "SINEC_Cf", "L1ME","L1MD", "L1MC", "L1MB", "L1MA", "L1_Carn", "L1_Canid", "L1_Canis", "L1_Cf", "L2", "MIR"),
                              repType = c("new_SINE", "new_SINE", "new_SINE", "new_SINE", "new_SINE", "old_LINE", "old_LINE", "old_LINE", "old_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "new_LINE", "ancient", "ancient")
                              )
  }
  
  return(repFamStruct)
}




