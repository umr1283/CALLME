cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

samplesFile <- paste(dir,"/sample_names.txt",sep="")
nFilesPerList <- as.numeric(nFilesPerList)

start <- 1
samples <- read.table(samplesFile,sep="\t",h=TRUE)
samples[,"Filename"] <- paste(dir,rep("/",nrow(samples)),samples[,"Filename"],sep="")
fileNum <- 1
lists <- character()
while(start<=nrow(samples)){
  stopp <- min(start+nFilesPerList,nrow(samples))
  write.table(samples[start:stopp,],paste(dir,"/list",fileNum,".txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
  start <- stopp+1
  fileNum <- fileNum + 1
  lists <- c(lists,paste(dir,"/list",fileNum-1,".txt",sep=""))
}
writeLines(lists,paste(dir,"/lists.txt",sep=""))
