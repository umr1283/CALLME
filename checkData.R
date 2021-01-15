cat(paste("Analysis started",date(),"\n"))
#dir <- "."
#setwd(dir)
cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

print(paste("data dir is",dir))

options(stringsAsFactors=FALSE)

setwd(dir)


indivs <- list.files("rawData")
indivs <- gsub(".txt","",indivs)

genocol <- 6
bafcol <- 5
log2ratiocol <- 4

IID <- character(0)
sdLRR <- numeric(0)
sdHetBAF <- numeric(0)
callRate <- numeric(0)
for(i in indivs){
  d <- read.table(paste("rawData/",i,".txt",sep=""),sep="\t",h=T)                 
  IID <- c(IID,i)
  sdLRR <- c(sdLRR,sd(d[,log2ratiocol],na.rm=TRUE))
  sdHetBAF <- c(sdHetBAF,sd(d[which(d[,genocol]=="AB"),bafcol],na.rm=TRUE))
  callRate <- c(callRate,1-length(which(d[,genocol]%in%c("NC","NoCall")))/nrow(d))
}

write.table(data.frame(IID,sdLRR,sdHetBAF,callRate),"data_check.txt",sep="\t",row.names=FALSE,quote=FALSE)

cat(paste("Analysis finished",date(),"\n"))

