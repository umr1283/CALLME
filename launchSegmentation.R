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

print(paste("data dir is",dir))

options(stringsAsFactors=FALSE)
library(gada)
Tcutoff <- as.numeric(Tcutoff)
minSegLength <- as.numeric(minSegLength)# ~500 KB
aAlpha <- as.numeric(aAlpha)
thresholdKB <- 2000
outdir <- paste(dir,"/results_T_",Tcutoff,"_aAlpha_",aAlpha,"_minSegLength_",minSegLength,sep="")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

setwd(dir)


indivs <- list.files("rawData")
indivs <- gsub(".txt","",indivs)

genocol <- 6
bafcol <- 5
log2ratiocol <- 4

if(!file.exists(paste(outdir,"/segments.txt",sep=""))){
  cat("#### DATA IMPORT ######\n\n")
  dat <- setupParGADA.B.deviation(NumCols=6,GenoCol=genocol, BAFcol=bafcol, log2ratioCol=log2ratiocol) 
  cat("#### SEGMENTATION ######\n\n")
  parSBL(dat, estim.sigma2=TRUE, aAlpha=aAlpha)
  cat("#### BACKWARD ELIMINATION ######\n\n")
  parBE.B.deviation(dat, T=Tcutoff, MinSegLen=minSegLength) 
  exportSegments2File(dat,file=paste(outdir,"/segments.txt",sep=""))
}
segments <- read.table(paste(outdir,"/segments.txt",sep=""),h=T,sep="\t")
segments <- segments[,c("IniProbe","EndProbe","LenProbe","chr","LRR","sample")]
segments$LenKB <- round((segments[,"EndProbe"]-segments[,"IniProbe"])/1000,3)
segments <- segments[which(segments[,"LenKB"]>thresholdKB),]
segments$SegID <- paste(segments[,"sample"],rep("_chr",nrow(segments)),segments[,"chr"],rep(":",nrow(segments)),segments[,"IniProbe"],rep("-",nrow(segments)),segments[,"EndProbe"],sep="")
Status <- character(nrow(segments))
names(Status) <- segments[,"SegID"]


cat("Assigning status to each event \n")
datas <- list()
k <- 1
for(ind in indivs){
  segInd <- segments[which(segments[,"sample"]==ind),]
  if(nrow(segInd)>0){
    datas[[ind]] <- read.table(paste("rawData/",ind,".txt",sep=""),h=T,sep="\t")
    for(i in 1:nrow(segInd)){
      dataSeg <- datas[[ind]][which(datas[[ind]][,"Chr"]==segments[i,"chr"] & datas[[ind]][,"Position"]>=segments[i,"IniProbe"] & datas[[ind]][,"Position"]<=segments[i,"EndProbe"]),]
      med <- median(dataSeg[,log2ratiocol],na.rm=TRUE)
      sd <- sd(dataSeg[,log2ratiocol],na.rm=TRUE)
      Status[segments[i,"SegID"]] <- "neutral"
      if(med>0.2*sd){
        Status[segments[i,"SegID"]] <- "gain"
      }
      else if(med < -0.2*sd){
        Status[segments[i,"SegID"]] <- "loss"
      }
    }
  }
}

segments$Status <- Status

segments <- segments[,c("IniProbe","EndProbe","LenProbe","chr","LRR","sample","LenKB","SegID","Status")]

merge2Segs <- function(segData){
  res <- data.frame(IniProbe=segData[1,"IniProbe"],
                    EndProbe=segData[nrow(segData),"EndProbe"],
                    LenProbe=sum(segData[,"LenProbe"]),
                    chr=segData [1,"chr"],
                    LRR=mean(segData[,"LRR"]),
                    sample=segData[1,"sample"],
                    LenKB=round((segData[nrow(segData),"EndProbe"]-segData[1,"IniProbe"])/1000,2),
                    SegID=paste(segData[1,"sample"],"_chr",segData[1,"chr"],":",segData[1,"IniProbe"],"-",segData[nrow(segData),"EndProbe"],sep="") ,
                    Status=segData[1,"Status"]
                    )
  #print(head(res))
  return(res)
}

cat("Merging adjacent events \n")
segMerged <- NULL
for(ind in indivs){
  segInd <- segments[which(segments[,"sample"]==ind),]
  if(nrow(segInd)>0){
    for(chr in unique(segInd[,"chr"])){
      segChr <- segments[which(segments[,"chr"]==chr),]
      if(nrow(segChr)>1){
        segChr <- segChr[order(segChr[,"IniProbe"]),]
        i=1
        while(i<nrow(segChr)){
          start <- i
          while(segChr[i+1,"IniProbe"]-segChr[i,"EndProbe"]<1000000 && segChr[i+1,"Status"]==segChr[i,"Status"] && i<nrow(segChr)){
            i <- i+1
          }
          if(i>start){
            segMerged <- rbind(segMerged,merge2Segs(segChr[start:i,]))
            i <- i+1
          }
          else{
            segMerged <- rbind(segMerged,segChr[i,])
            i <- i+1
          }
        }
      }
      else{
        segMerged <- rbind(segMerged,segChr)
      }
    }
  }
}

if(is.null(segMerged)){
  segMerged <- segments
}

cat("Calculating number of BAF bands and percentage of mosaic cells \n")

pctCells <- rep(NA,nrow(segMerged))
names(pctCells) <- segMerged[,"SegID"]
mu1<- numeric(nrow(segMerged))
names(mu1) <- segMerged[,"SegID"]
mu2<- numeric(nrow(segMerged))
names(mu2) <- segMerged[,"SegID"]
isMOSAIC<- rep("No",nrow(segMerged))
names(isMOSAIC) <- segMerged[,"SegID"]
pctProbesG1 <- rep(NA,nrow(segMerged))
names(pctProbesG1) <- segMerged[,"SegID"]
pctProbesG2 <- rep(NA,nrow(segMerged))
names(pctProbesG2) <- segMerged[,"SegID"]

G <- 4
for(ind in indivs){
  segInd <- segMerged[which(segMerged[,"sample"]==ind),]
  if(nrow(segInd)>0){
    for(iEvent in 1:nrow(segInd)){
       select <- (datas[[ind]][,"Position"]>=segInd[iEvent,"IniProbe"] & datas[[ind]][,"Position"]<=segInd[iEvent,"EndProbe"] & datas[[ind]][,"Chr"]==segInd[iEvent,"chr"])
       tmp    <-  datas[[ind]][select,]
       x <- tmp[,bafcol]
       x <- x[!is.na(x)]
       #clus <- Mclust(x,G=G)
       #h <- function(k){ ain = kmeans(x,k,nstart=10); ain$betweenss/ain$totss }
       #hx <- sapply(Gset,h)
       #gKmeans <- Gset[which.max(hx)]
       #nBAFbands1[segInd[iEvent,"SegID"]] <- clus$G
       #nBAFbands2[segInd[iEvent,"SegID"]] <- gKmeans
       km <- kmeans(x,G,nstart=10)
       m <- sort(km$centers[,1])
       #m2 <- sort(kmeans(x,gKmeans,nstart=10)$centers[,1])
       m1 <-  m[2]
       m2 <- m[3]
       pctProbesG1[segInd[iEvent,"SegID"]] <- table(km$cluster)[names(m1)]/nrow(tmp) 
       pctProbesG2[segInd[iEvent,"SegID"]] <- table(km$cluster)[names(m2)]/nrow(tmp) 
       mu1[segInd[iEvent,"SegID"]] <- m1
       mu2[segInd[iEvent,"SegID"]] <- m2
       if(segInd[iEvent,"Status"]=="gain"){
          pctCells[segInd[iEvent,"SegID"]] <- 2*(m2-m1)/(1-m2+m1)  
       }
       else if(segInd[iEvent,"Status"]=="loss"){
          pctCells[segInd[iEvent,"SegID"]] <- 2*(m2-m1)/(1+m2-m1)  
       }
       else{
          pctCells[segInd[iEvent,"SegID"]] <- m2-m1             
       }
       dist1 <- 0.5 - m1
       dist2 <- m2 -0.5
       ratio <- min(c(dist1,dist2))/ max(c(dist1,dist2))
       if(m1<0.5 & m2>0.5 & pctCells[segInd[iEvent,"SegID"]]>0.05 &  pctCells[segInd[iEvent,"SegID"]]<0.95 & ratio > 0.5 & pctProbesG1[segInd[iEvent,"SegID"]]>0.05 & pctProbesG2[segInd[iEvent,"SegID"]]>0.05){
          isMOSAIC[segInd[iEvent,"SegID"]] <- "Yes"
       }
    }
  }
}
Datafile <- paste(rep(dir,nrow(segMerged)),rep("/",nrow(segMerged)),rep("/rawData/",nrow(segMerged)),segMerged[,"sample"],rep(".txt",nrow(segMerged)),sep="")


segMerged <- cbind(segMerged,pctCells,mu1,mu2,pctProbesG1,pctProbesG2,isMOSAIC,Datafile)

write.table(segMerged,paste(outdir,"/segments.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
if(file.exists("SBL")){
  system(paste("mv SBL",outdir))
}
cat(paste("Analysis finished",date(),"\n"))

stop()
cat("#### PLOTTING RESULTS  ######\n\n")
labels <- c("gain","loss","neutral")
col1 <- rgb(0,0,0,alpha=0.5)
cols <- c(rgb(1,0,0,alpha=0.25),rgb(0,0,1,alpha=0.25),rgb(0,1,0,alpha=0.25))
names(cols) <- labels
for(ind in indivs){
  cat(paste("  - Creating plots for",ind,"(",k,"/",length(indivs),")\n"))
  segInd <- segMerged[which(segments[,"sample"]==ind),]
  if(nrow(segInd)>0){
    dataInd <- datas[[ind]]
      for(chr in 1:22){
      segChr <- segInd[which(segInd[,"chr"]==chr),]
      if(nrow(segChr)>0){
        dataChr <- dataInd[which(dataInd[,"Chr"]==chr),]
        for(i in 1:nrow(segChr)){
           suffix <- ""
           if(segChr[i,"isMOSAIC"]=="Yes"){
             suffix="_mosaic"
           }
           png(paste(outdir,"/chr",chr,"_",segChr[i,"IniProbe"],"-",segChr[i,"EndProbe"],suffix,".png",sep=""))
           xlim1 <- max(min(segChr[i,"IniProbe"]-0.5*segChr[i,"LenKB"]*1000,segChr[i,"IniProbe"]-2000000),0)
           xlim2 <- min(max(segChr[i,"EndProbe"]+0.5*segChr[i,"LenKB"]*1000,segChr[i,"EndProbe"]+2000000),max(dataChr[,"Position"],na.rm=TRUE))
           dataSeg <- dataChr[which(dataChr[,"Position"]>=segChr[i,"IniProbe"] & dataChr[,"Position"]<=segChr[i,"EndProbe"]),] 
           xlim <- c(xlim1,xlim2)
           toPlot <- dataChr[which(dataChr[,"Position"]>=xlim[1] & dataChr[,"Position"]<=xlim[2]),]                                   
           layout(matrix(c(1,0,2),nrow=3),heights=c(0.45,0.05,0.45))
           plot(toPlot[,3],toPlot[,log2ratiocol],col=col1,ylim=c(-2,2),xlim=xlim,pch=20,cex=0.5,ylab="LRR",xlab="Position")
           abline(h=0,lwd=0.3,col="darkgrey")
           rect(segChr[i,"IniProbe"],-2,segChr[i,"EndProbe"],2,col=cols[as.character(segChr[i,"Status"])],pch=20,cex=0.5)
           title(main=paste("Indiv ",ind," (chr ",chr,")\n",segChr[i,"Status"],", ",round(segChr[i,"pctCells"],2)," affected, mosaic : ",segChr[i,"isMOSAIC"],"\n(",round(segChr[i,"LenKB"])," KB)",sep=""))
           plot(toPlot[,3],toPlot[,bafcol],col=col1,ylim=c(0,1),xlim=xlim,pch=20,cex=0.5,ylab="BAF",xlab="Position")
           abline(h=0.5,lwd=0.3,col="darkgrey")
           abline(h=1/3,lwd=0.3,col="darkgrey")
           abline(h=2/3,lwd=0.3,col="darkgrey")
           abline(h=c(segChr[i,"mu1"],segChr[i,"mu2"]),col="red",lty=2)
	   rect(segChr[i,"IniProbe"],0,segChr[i,"EndProbe"],1,col=cols[as.character(segChr[i,"Status"])],pch=20,cex=0.5)
  	   title(main=paste("Indiv ",ind," (chr ",chr,")\n",segChr[i,"Status"],", ",round(segChr[i,"pctCells"],2)," affected, mosaic : ",segChr[i,"isMOSAIC"],"\n(",round(segChr[i,"LenKB"])," KB)",sep=""))
           dev.off()
        }
      }
    }
  }
  k <- k+1

}

if(file.exists("SBL")){
  system(paste("mv SBL",outdir))
}
cat(paste("Analysis finished",date(),"\n"))

