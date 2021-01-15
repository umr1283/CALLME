cat(paste("Analysis started",date(),"\n"))

cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

options(stringsAsFactors=FALSE)

indivs <- readLines(indivsFile)

checkALL <- NULL
for(i in indivs){
    check <- read.table(paste(dir,"/",i,"/data_check.txt",sep=""),h=T,sep="\t")
    checkALL <- rbind(checkALL,check)
}

write.table(checkALL,paste(dir,"/data_check_all.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)


cat(paste("Analysis finished",date(),"\n"))
