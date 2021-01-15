cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

nFilesPerList <- as.numeric(nFilesPerList)

start <- 1
files <- list.files(dir,full.names=TRUE)
files <- files[grep("PennCNV",files)]
fileNum <- 1
dir.create(paste(dir,"/GC_corrected/logs",sep=""),recursive=TRUE,showWarnings=FALSE)

while(start<=length(files)){
  stopp <- min(start+nFilesPerList,length(files))
  writeLines(files[start:min(start+nFilesPerList,length(files))],paste(dir,"/list",fileNum,".txt",sep=""))
  system(paste("nohup perl ",libDir,"/genomic_wave.pl --adjust --gcmodelfile=",libDir,"/hhall.hg18.gcmodel --listfile=",dir,"/list",fileNum,".txt --prefix=",dir,"/GC_corrected >& ",dir,"/GC_corrected/logs/correction_",fileNum,".out &",sep=""))
  start <- stopp+1
  fileNum <- fileNum + 1
}
