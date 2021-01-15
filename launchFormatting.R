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
files <- list.files(dir2,full.names=TRUE)
files <- files[grep("_extracted",files)]
#filesToExclude <- readLines(paste(dir2,"/toExclude.txt",sep=""))
#filesToExclude <- paste(dir2,"/",filesToExclude,"_extracted.txt",sep="")
#files <- files[which(!files%in%filesToExclude)]
              
fileNum <- 1
dir.create(paste(outdir,"/logs",sep=""),recursive=TRUE,showWarnings=FALSE)
while(start<=length(files)){
  stopp <- min(start+nFilesPerList,length(files))
  writeLines(files[start:min(start+nFilesPerList,length(files))],paste(dir2,"/list",fileNum,".txt",sep=""))
  system(paste("nohup perl /ep10/disks/SANA8/boris/home/scripts/MAD/formatData3.pl ",dir1," 3 0 ",dir2,"/list",fileNum,".txt ",outdir," >& ",outdir,"/logs/formatData",fileNum,".out &",sep=""))
  start <- stopp+1
  fileNum <- fileNum + 1
}
