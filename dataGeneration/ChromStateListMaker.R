#### so if we get histone stuff as a rep file 

### we could probbly read in each red as though it were a family 

### set it up like the repdb?  or something differnt 

### chose cell line and extract regions 
### alternativly i could load a bunch of objects 

Lines <- list.files("~/Desktop/Domain_manuscript/Data/hmmCombined/")
Lines <- gsub(pattern = "_combined.txt", replacement = "",x = Lines)

for(i in 1:length(Lines)){
  cellLine = Lines[i]
  cl <- read.table(file = paste("~/Desktop/Domain_manuscript/Data/hmmCombined/", cellLine, "_combined.txt", sep = ""), 
                   header = TRUE, comment.char = "",
                   colClasses = c("integer", "character", "integer", "integer", "factor", "integer", "character", "integer", "integer", "character"))
  # get the cell line and split the file into regions 
  #assign it at the end and save all the lists in an R object
  clList <- split(x = cl,f = cl$name)
  assign(x = cellLine , clList)
}

save(list = Lines,file = "~/Desktop/Domain_manuscript/R_objects/chromStateCombined")











