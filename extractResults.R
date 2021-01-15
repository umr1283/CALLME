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

T <- as.numeric(T)
minSegLength <- as.numeric(minSegLength)
aAlpha <- as.numeric(aAlpha)
outdir <- paste(dir,"/results_T_",T,"_aAlpha_",aAlpha,"_minSegLength_",minSegLength,sep="")
indivs <- readLines(indivsFile)
sdBAFthreshold <- as.numeric(sdBAFthreshold)
sdLRRthreshold <- as.numeric(sdLRRthreshold)
toExclude <- character(0)
if(file.exists(paste(dir,"/data_check_all.txt",sep=""))){
  check <- read.table(paste(dir,"/data_check_all.txt",sep=""),h=T,sep="\t")
  toExclude <- check[which(check$sdHetBAF>sdBAFthreshold | check$sdLRR>sdLRRthreshold ),"IID"]
  cat(paste(length(toExclude),"samples excluded due to unsufficient data quality\n"))
}


thresholdKB <- 2000

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)


cat("Reading all segments\n")
segmentsall <- NULL
for(i in indivs){
  if(!i%in%toExclude){
    if(file.exists(paste(dir,"/",i,"/results_T_",T,"_aAlpha_",aAlpha,"_minSegLength_",minSegLength,"/segments.txt",sep=""))){
      seg <- read.table(paste(dir,"/",i,"/results_T_",T,"_aAlpha_",aAlpha,"_minSegLength_",minSegLength,"/segments.txt",sep=""),h=T,sep="\t")
      segmentsall <- rbind(segmentsall,seg)
    }
    else{
      cat(paste("WARNING : file",paste(dir,"/",i,"/results_T_",T,"_aAlpha_",aAlpha,"_minSegLength_",minSegLength,"/segments.txt",sep=""),"does not exist"))
    }
  }
}

write.table(segmentsall,paste(outdir,"/segments_all.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
cat(paste(" => ",nrow(segmentsall)," found for ",length(indivs),"individuals (see in ",outdir,"/segments_all.txt)\n\n",sep=""))

possiblyMosaics <- segmentsall[which(segmentsall[,"isMOSAIC"]=="Yes" & segmentsall[,"LenKB"]>thresholdKB),]



cat("Plotting mosaics\n")
indchr <- paste(possiblyMosaics[,"sample"],possiblyMosaics[,"chr"],sep="_")
possiblyMosaics$IndChr <- indchr
segms <- possiblyMosaics
genocol <- 6
bafcol <- 5
log2ratiocol <- 4

labels <- c("gain","loss","neutral")
col1 <- rgb(0,0,0,alpha=0.5)
cols <- c(rgb(1,0,0,alpha=0.25),rgb(0,0,1,alpha=0.25),rgb(0,1,0,alpha=0.25))
names(cols) <- labels

k <- 1
pdf(paste(outdir,"/possible_mosaics.pdf",sep=""))
for(ind in unique(possiblyMosaics[,"sample"])){
    cat(paste("Processing sample",ind,"(",k,"/",length(unique(possiblyMosaics[,"sample"])),")\n"))
    dataInd <- read.table(paste(dir,"/",ind,"/rawData/",ind,".txt",sep=""),h=T,sep="\t")
    segs <- segms[which(segms[,"sample"]==ind),]
    for(chr in unique(segs[,"chr"])){
           forRect <- possiblyMosaics[which(possiblyMosaics[,"IndChr"]==paste(ind,chr,sep="_")),]
           dataChr <- dataInd[which(dataInd[,"Chr"]==chr),]
           xlim <- c(min(dataChr[,"Position"]),max(dataChr[,"Position"]))
           toPlot <- dataChr                                   
           layout(matrix(c(1,0,2),nrow=3),heights=c(0.45,0.05,0.45))
           plot(toPlot[,3],toPlot[,log2ratiocol],col=col1,ylim=c(-2,2),xlim=xlim,pch=20,cex=0.5,ylab="LRR",xlab="Position")
           abline(h=0,lwd=0.3,col="darkgrey")
           rect(forRect[,"IniProbe"],-2,forRect[,"EndProbe"],2,col=cols[as.character(forRect[,"Status"])],pch=20,cex=0.5)
  	   title(main=paste("Indiv ",ind," (chr ",chr,")\n",toString(forRect[,"LenKB"])," KB || ",toString(forRect[,"Status"])," || ",toString(round(forRect[,"pctCells"]*100))," % affected\n",sep=""))
           plot(toPlot[,3],toPlot[,bafcol],col=col1,ylim=c(0,1),xlim=xlim,pch=20,cex=0.5,ylab="BAF",xlab="Position")
           abline(h=0.5,lwd=0.3,col="darkgrey")
           abline(h=1/3,lwd=0.3,col="darkgrey")
           abline(h=2/3,lwd=0.3,col="darkgrey")
           abline(h=c(segs[i,"mu1"],segs[i,"mu2"]),col="red",lty=2)
	   rect(forRect[,"IniProbe"],0,forRect[,"EndProbe"],1,col=cols[as.character(forRect[,"Status"])],pch=20,cex=0.5)
  	   title(main=paste("Indiv ",ind," (chr ",chr,")\n",toString(forRect[,"LenKB"])," KB || ",toString(forRect[,"Status"])," || ",toString(round(forRect[,"pctCells"]*100))," % affected\n",sep=""))

     }
     k <- k+1
  }

dev.off()

cat(paste("Analysis finished",date(),"\n"))

