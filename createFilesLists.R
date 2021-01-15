dir <- "/ep10/disks/SANA8/boris/MOSAICS/METABO_DESIR/normalization/normalized"
nFilesPerList <- 100

start <- 1
files <- list.files(dir,full.names=TRUE)
files <- files[grep("PennCNV",files)]
fileNum <- 1
while(start<=length(files)){
  stopp <- min(start+nFilesPerList,length(files))
  writeLines(files[start:min(start+nFilesPerList,length(files))],paste(dir,"/list",fileNum,".txt",sep=""))
  system(paste("perl /ep10/disks/SANA8/boris/home/scripts/MAD/genomic_wave.pl --adjust --gcmodelfile=/ep10/disks/SANA8/boris/home/bin/penncnv/lib/hhall.hg18.gcmodel --listfile=",dir,"/list",fileNum,".txt --prefix=",dir,"/GC_corrected >list",fileNum,".out &",sep=""))
  start <- stopp+1
  fileNum <- fileNum + 1
}
