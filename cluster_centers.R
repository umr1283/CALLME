cat(paste("Analysis started",date(),"\n\n"))
options(stringsAsFactors=TRUE)
cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

### PARAMETERS
# samples.file
# input.path
# geno.path
# out.file



### Read the file with assays to process ####
# file name should be 'sample_names.txt' 
# format: Assay
assaysToRow<-read.table(samples.file,sep="\t",header=TRUE)
print(assaysToRow)
#print(head(assaysToRow))
nbrAssays<-length(assaysToRow$Assay)
###


### Read the appropiate tQN cluster file
# format:
#reporterId	AA_T_Mean	AA_T_Dev	AB_T_Mean	AB_T_Dev	BB_T_Mean	BB_T_Dev	AA_R_Mean	AA_R_Dev	AB_R_Mean	AB_R_Dev	BB_R_Mean	BB_R_Dev


### Require the R package Limma from www.bioconductor.org
library(limma)
###

cat("Performing normalization\n")
## perform tQN normalization ###
for(r in 1:nbrAssays){
	#foreach master assay, calculate, write and plot for the 4 different platforms.
	sampleName<-assaysToRow$Assay[r]
	cat(paste("  - sample",sampleName," (",r,"/",nbrAssays,")\n"))
	
	#### Read sample data ####
	my.file<-paste(input.path,sampleName,"_extracted.txt",sep="")
	baf.data<-read.delim(my.file,header=TRUE,na.strings=c(NA,NaN))
	baf.data$Name<-as.character(baf.data$Name)
        #format (tab separated):
	#Name	Chr	Position	X	Y
	my.file<-paste(geno.path,sampleName,"_extracted.txt",sep="")
	geno.data<-read.delim(my.file,header=TRUE,na.strings=c(NA,NaN))
	geno.data$Name<-as.character(baf.data$Name)
        
        baf.data <- merge(baf.data,geno.data,by.x="Name",by.y="Name")
        print(head(baf.data))
	###

	### Check for presence of CNV probes, as these should not be tQN normalized ###
	cnv_present<-FALSE
	uu<-grep("cnv",baf.data$Name)
	cnv.data<-data.frame()
	if(length(uu)>0){
		cnv_present<-TRUE
		cnv.data<-baf.data[uu,]
		baf.data<-baf.data[-uu,]
	}

        
	#### Quantile normalization ####
	cat("      * Quantile normalization \n")
	AA<-normalizeQuantiles(cbind(baf.data$X,baf.data$Y))
	####
	
	#### Collect R ####
	R.tQN<-AA[,1]+AA[,2] # Rvalue = int Y + int X
	R.cnv<-c()
	if(cnv_present){
		R.cnv<-cnv.data$X+cnv.data$Y
	}
	####
	
	
	#### Thresholding ####
	QN.effect.X<-AA[,1]/baf.data$X
	x.threshold<-1.5
	aff.x<-which(QN.effect.X>x.threshold)
		
	QN.effect.Y<-AA[,2]/baf.data$Y
	y.threshold<-1.5
	aff.y<-which(QN.effect.Y>y.threshold)

	if(length(aff.x)>0){
		AA[aff.x,1]<-x.threshold*baf.data$X[aff.x]
	}
	if(length(aff.y)>0){
		AA[aff.y,2]<-y.threshold*baf.data$Y[aff.y]
	}
	####
	
	
	#### Calculate theta ###
	theta.tQN<-2/pi*atan(AA[,2]/AA[,1])
	theta.cnv<-c()
	if(cnv_present){
		theta.cnv<-2/pi*atan(cnv.data$Y/cnv.data$X)
	}
	####
	
	if(r==1){
	    nAA <- numeric(nrow(baf.data))
	    nAA[which(baf.data[,"GType"]=="AA")] <- 1 
	    TAA <- numeric(nrow(baf.data))
	    TAA[which(baf.data[,"GType"]=="AA")] <- theta.tQN[which(baf.data[,"GType"]=="AA")]
	    RAA <- numeric(nrow(baf.data))
	    RAA[which(baf.data[,"GType"]=="AA")] <- R.tQN[which(baf.data[,"GType"]=="AA")]
            nAB <- numeric(nrow(baf.data))
	    nAB[which(baf.data[,"GType"]=="AB")] <- 1 
	    TAB <- numeric(nrow(baf.data))
	    TAB[which(baf.data[,"GType"]=="AB")] <- theta.tQN[which(baf.data[,"GType"]=="AB")]
	    RAB <- numeric(nrow(baf.data))
	    RAB[which(baf.data[,"GType"]=="AB")] <- R.tQN[which(baf.data[,"GType"]=="AB")]
	    nBB <- numeric(nrow(baf.data))
	    nBB[which(baf.data[,"GType"]=="BB")] <- 1 
	    TBB <- numeric(nrow(baf.data))
	    TBB[which(baf.data[,"GType"]=="BB")] <- theta.tQN[which(baf.data[,"GType"]=="BB")]
	    RBB <- numeric(nrow(baf.data))
	    RBB[which(baf.data[,"GType"]=="BB")] <- R.tQN[which(baf.data[,"GType"]=="BB")]

      }
      else{
	    nAA[which(baf.data[,"GType"]=="AA")] <- nAA[which(baf.data[,"GType"]=="AA")] + 1 
	    TAA[which(baf.data[,"GType"]=="AA")] <- TAA[which(baf.data[,"GType"]=="AA")] + theta.tQN[which(baf.data[,"GType"]=="AA")]
	    RAA[which(baf.data[,"GType"]=="AA")] <- RAA[which(baf.data[,"GType"]=="AA")] + R.tQN[which(baf.data[,"GType"]=="AA")]
            nAB[which(baf.data[,"GType"]=="AB")] <- nAB[which(baf.data[,"GType"]=="AB")] + 1 
	    TAB[which(baf.data[,"GType"]=="AB")] <- TAB[which(baf.data[,"GType"]=="AB")] + theta.tQN[which(baf.data[,"GType"]=="AB")]
	    RAB[which(baf.data[,"GType"]=="AB")] <- RAB[which(baf.data[,"GType"]=="AB")] + R.tQN[which(baf.data[,"GType"]=="AB")]
	    nBB[which(baf.data[,"GType"]=="BB")] <- nBB[which(baf.data[,"GType"]=="BB")] + 1 
	    TBB[which(baf.data[,"GType"]=="BB")] <- TBB[which(baf.data[,"GType"]=="BB")] + theta.tQN[which(baf.data[,"GType"]=="BB")]
	    RBB[which(baf.data[,"GType"]=="BB")] <- RBB[which(baf.data[,"GType"]=="BB")] + R.tQN[which(baf.data[,"GType"]=="BB")]

	}
}

meanTAA <- TAA/nAA
meanTAA[is.infinite(TAA)] <- NA
meanTAB <- TAB/nAB
meanTAB[is.infinite(TAB)] <- NA
meanTBB <- TBB/nBB
meanTBB[is.infinite(TBB)] <- NA
meanRAA <- RAA/nAA
meanRAA[is.infinite(RAA)] <- NA
meanRAB <- RAB/nAB
meanRAB[is.infinite(RAB)] <- NA
meanRBB <- RBB/nBB
meanRBB[is.infinite(RBB)] <- NA

res <- data.frame(reporterId=baf.data[,"Name"],
			AA_T_Mean=meanTAA,
			AB_T_Mean=meanTAB,
			BB_T_Mean=meanTBB,
			AA_R_Mean=meanRAA,
			AB_R_Mean=meanRAB,
			BB_R_Mean=meanRBB)
write.table(res,out.file,sep="\t",row.names=FALSE,quote=FALSE)



cat(paste("Analysis finished",date(),"\n\n"))

#end each master assay	
