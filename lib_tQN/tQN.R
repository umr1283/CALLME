# Copyright (C) 2008 Johan Staaf
# 
# This file is part of tQN,
# http://baseplugins.thep.lu.se/wiki/se.lu.onk.IlluminaSNPNormalization
# 
# tQN is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
# 
# tQN is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Class Discoverer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307 USA


### This scripts performs tQN normalization of Illumina Infinium II data
#   Files to be processed should be listed in a txt file called sample_names.txt
#   Format of sample_names.txt should be a single column of sample names with the header:
#   Assay
#   For each sample a separate file should exist. The file name should be the sample name with
#   suffix  _extracted.txt
#   sample_names.txt and individual sample files should be in the same directory.

#   During normalization, only SNPs present in the tQN normalized cluster files will be propageted. Remaining SNPs
#   will be lost as no cluster assignment is obtainable.


cat(paste("Analysis started",date(),"\n\n"))
### File paths and other parameters
#cluster.path<-"/ep10/disks/SANA8/boris/home/bin/tQN-1.1.2/lib/"
#source("tQN_parameters.txt")
cmdline <- commandArgs()
for (e in cmdline[-(1:2)]){
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2])){
    assign(ta[[1]][1],as.character(ta[[1]][2]))
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

###


### Read the file with assays to process ####
# file name should be 'sample_names.txt' 
# format: Assay
assaysToRow<-read.table(samples.file,sep="\t",header=TRUE)
print(head(assaysToRow))
#print(head(assaysToRow))
nbrAssays<-length(assaysToRow$Assay)
###


### Read the appropiate tQN cluster file
# format:
#reporterId	AA_T_Mean	AA_T_Dev	AB_T_Mean	AB_T_Dev	BB_T_Mean	BB_T_Dev	AA_R_Mean	AA_R_Dev	AB_R_Mean	AB_R_Dev	BB_R_Mean	BB_R_Dev
my.cluster.file<-bead.platform
my.cluster.file<-paste(cluster.path,my.cluster.file,sep="")
cat("Reading cluster file\n")
clusterFile.new<-read.delim(my.cluster.file,header=TRUE)
clusterFile.new$reporterId<-as.character(clusterFile.new$reporterId)
###


### BAF Function ###
baf.function<-function(Theta,tAA,tAB,tBB){
	BAF<-Theta
	tmp.tBB<-median(tBB,na.rm=TRUE)
	tmp.tAA<-median(tAA,na.rm=TRUE)
        for(t in 1:length(Theta)){
		if(is.na(Theta[t])){
			BAF[t]<-NaN
			next
		}
				
		e.tAA<-TRUE
		e.tAB<-TRUE
		e.tBB<-TRUE
		
		if(is.na(tAA[t])){
			e.tAA<-FALSE
		}
		if(is.na(tAB[t])){
			e.tAB<-FALSE
		}
		if(is.na(tBB[t])){
			e.tBB<-FALSE
		}
		
		# 0: Test for inconsistencies between tAA, tAB, tBB
		if( e.tAA & e.tAB){
			if(tAA[t]>tAB[t]){
				# tAA has a higher value than tAB which is inconsistent! BLANK SNP
				BAF[t]<-NaN
				next
			}
		}
		if( e.tAA & e.tBB){
			if(tAA[t]>tBB[t]){
				# tAA has a higher value than tBB which is inconsistent! BLANK SNP
				BAF[t]<-NaN
				next
			}
		}
		if( e.tAB & e.tBB){
			if(tAB[t]>tBB[t]){
				# tAB has a higher value than tBB which is inconsistent! BLANK SNP
				BAF[t]<-NaN
				next
			}
		}
		
		
		# 1: Triple blank SNP
		if(!(e.tAA) & !(e.tAB) & !(e.tBB)){
			# no GType clusters exist for this SNP. BLANK IT!
			BAF[t]<-NaN
		}
		
		# 2: Blank for AB, AA, while positive for BB
		if(!(e.tAA) & !(e.tAB) & (e.tBB)){
			# no AA & AB cluster exist for this SNP, while BB exists. Set it to BB
			if(Theta[t]>=tBB[t]){
				BAF[t]<-1
			}else{
				BAF[t]<-NaN
			}
		}
		
		# 3: Blank for AB, BB, while positive for AA
		if((e.tAA) & !(e.tAB) & !(e.tBB)){
			# no BB & AB cluster exist for this SNP, while AA exists. Set it to AA
			if(Theta[t]<= tAA[t]){
				BAF[t]<-0
			}else{
				BAF[t]<-NaN
			}
		}

		# 4: Blank for AB while positive for AA & BB
		if((e.tAA) & !(e.tAB) & (e.tBB)){
			# no AB cluster exist for this SNP, while AA & BB exists. Set it to the closest of AA or BB
			min.index<-which.min(c(abs(tAA[t]-Theta[t]),abs(tBB[t]-Theta[t]) ) )
			if(min.index==1){
				# closest to tAA
				if(Theta[t]< tAA[t]){
					BAF[t]<-0
				}else{
					BAF[t]<-NaN
				}
			}else{
				if(Theta[t]>=tBB[t]){
					BAF[t]<-1
				}else{
					BAF[t]<-NaN
				}
			}
		}
		
		# 5: Blank for AA while positive for AB & BB
		if( !(e.tAA) & (e.tAB) & (e.tBB)){
			if(Theta[t]>=tBB[t]){
				BAF[t]<-1
			}else{
				#two options!
				# 1: SNP is "correctly between" tAB and tBB
				# 2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be between tAA and tAB where it normally should not NOT BE.
				if(Theta[t]>=tAB[t]){
					#option 1! the expected
					#interpolate as SNP is expected to be between tAB and tBB
					BAF[t]<- 0.5+0.5*(Theta[t]-tAB[t])/(tBB[t]-tAB[t])
				}else{
					#option 2, NOT expected
					
					if(Theta[t]<tmp.tAA){
						#below tmp.tAA, set to 0
						BAF[t]<-0
					}else{
						#between tmp.tAA and tAB, interpolate
						BAF[t]<- 0.5*(Theta[t]-tmp.tAA)/(tAB[t]-tmp.tAA) 
					}
					
				}
			}
		}
		
		# 6: Blank for BB while positive for AA & AB
		if( (e.tAA) & (e.tAB) & !(e.tBB)){
			if(Theta[t]< tAA[t]){
				BAF[t]<-0
			}else{
				#two options!
				# 1: SNP is "correctly between" tAA and tAB
				# 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be between tAB and tBB where it normally should not NOT BE.
				if(Theta[t]<=tAB[t]){
					#option 1! the expected
					#interpolate as SNP is expected to be between tAB and tBB
					BAF[t]<- 0.5*(Theta[t]-tAA[t])/(tAB[t]-tAA[t])
				}else{
					#option 2, NOT expected
					
					if(Theta[t]>tmp.tBB){
						#above tmp.tBB, set to 1
						BAF[t]<-1
					}else{
						#between tAB and tmp.tBB, interpolate
						BAF[t]<- 0.5+0.5*(Theta[t]-tAB[t])/(tmp.tBB[t]-tAB[t])
					}
				}
			}
		}
		
		# 7: positive for AA & BB & AB
		if((e.tAA) & (e.tAB) & (e.tBB)){
			# AA & BB & AB exists. DO Illumina style calculation
			if(Theta[t]<tAB[t]){
				#lower part of BAF scale (0-0.5 typically)
				if(Theta[t]< tAA[t]){
					BAF[t]<-0
				}else{
					#interpolate as SNP is between tAA and tAB
					BAF[t]<- 0.5*(Theta[t]-tAA[t])/(tAB[t]-tAA[t])
				}
			}else{
				#upper part of BAF scale (0.5-1 typically)
				if(Theta[t]>=tBB[t]){
					BAF[t]<-1
				}else{
					#interpolate as SNP is between tAB and tBB
					BAF[t]<- 0.5+0.5*(Theta[t]-tAB[t])/(tBB[t]-tAB[t])
				}
			}
		}
	} #end for
	rm(tmp.tBB)
	rm(tmp.tAA)
	
	uu<-which(BAF>1)
	if(length(uu)>0){
		BAF[uu]<-1
	}
	uu<-which(BAF<0)
	if(length(uu)>0){
		BAF[uu]<-0
	}
	rm(uu)
	return(BAF)
}
###


### LogR ratio Function ###
logR.function<-function(R,Theta,tAA,tAB,tBB,rAA,rAB,rBB){
	#eR=rAA+(theta-tAA)*(rAB-rAA)/(tAB-tAA) if tAA<theta<tAB
	#eR=rAB+(theta-tAB)*(rBB-rAB)/(tBB-tAB) if tAB<theta<tBB
	#logR = log2(R/eR)
	logR<-R
	tmp.tBB<-median(tBB,na.rm=TRUE)
	tmp.tAA<-median(tAA,na.rm=TRUE) 
	for(t in 1:length(Theta)){
		if(is.na(Theta[t])){
			logR[t]<-NaN
			next
		}
		e.tAA<-TRUE
		e.tAB<-TRUE
		e.tBB<-TRUE
		
		if(is.na(tAA[t])){
			e.tAA<-FALSE
		}
		if(is.na(tAB[t])){
			e.tAB<-FALSE
		}
		if(is.na(tBB[t])){
			e.tBB<-FALSE
		}
		
		# 0: Test for inconsistencies between tAA, tAB, tBB
		if( e.tAA & e.tAB){
			if(tAA[t]>tAB[t]){
				# tAA has a higher value than tAB which is inconsistent! BLANK SNP
				logR[t]<-NaN
				next
			}
		}
		if( e.tAA & e.tBB){
			if(tAA[t]>tBB[t]){
				# tAA has a higher value than tBB which is inconsistent! BLANK SNP
				logR[t]<-NaN
				next
			}
		}
		if( e.tAB & e.tBB){
			if(tAB[t]>tBB[t]){
				# tAB has a higher value than tBB which is inconsistent! BLANK SNP
				logR[t]<-NaN
				next
			}
		}
		
		
		# 1: Triple blank SNP
		if(!(e.tAA) & !(e.tAB) & !(e.tBB)){
			# no GT cluster exist for this SNP. BLANK IT!
			logR[t]<-NaN
		}
		
		# 2: Blank for AB, AA, while positive for BB
		if(!(e.tAA) & !(e.tAB) & (e.tBB)){
			# no AA & AB cluster exist for this SNP, while BB exists. Set it to BB if Theta >= tBB
			#eR<-rBB[t]+(tBB[t]-Theta[t])*(rBB[t])/(tBB[t])
			if(Theta[t]>=tBB[t]){
				eR<-rBB[t]
				logR[t] <- log2(R[t]/eR)
			}else{
				logR[t]<-NaN
			}
		}
		
		# 3: Blank for AB, BB, while positive for AA
		if((e.tAA) & !(e.tAB) & !(e.tBB)){
			# no BB & AB cluster exist for this SNP, while AA exists. Set it to AA if Theta <= tAA
			if(Theta[t]<=tAA[t]){
				eR<-rAA[t]
				#eR<-rAA[t]+(Theta[t]-tAA[t])*(rAA[t])/(tAA[t])
				logR[t] <- log2(R[t]/eR)
			}else{
				logR[t]<-NaN
			}
		}

		# 4: Blank for AB while positive for AA & BB
		if((e.tAA) & !(e.tAB) & (e.tBB)){
			# no AB cluster exist for this SNP, while AA & BB exists. Set it to the closest of AA or BB
			min.index<-which.min(c(abs(tAA[t]-Theta[t]),abs(tBB[t]-Theta[t]) ) )
			if(min.index==1){
				# closest to tAA
				if(Theta[t]< tAA[t]){
					eR<-rAA[t]
					logR[t] <- log2(R[t]/eR)
				}else{
					logR[t]<-NaN
				}
			}else{
				#closest to tBB
				if(Theta[t]>=tBB[t]){
					eR<-rBB[t]
					logR[t] <- log2(R[t]/eR)
				}else{
					logR[t]<-NaN
				}
				
				
			}
		}
		
		# 5: Blank for AA while positive for AB & BB
		if( !(e.tAA) & (e.tAB) & (e.tBB)){
			
			#two options!
			# 1: SNP is "correctly between" tAB and tBB
			# 2: Heterozygous SNP is subjected to deletion or UPD of allele B making it unexectedly to be between tAA and tAB where it normally should not NOT BE.
			if(Theta[t]>=tAB[t]){
				#option 1! the expected
				#interpolate as SNP is expected to be between tAB and tBB
				eR<-rAB[t]+(Theta[t]-tAB[t])*(rBB[t]-rAB[t])/(tBB[t]-tAB[t])
				logR[t] <- log2(R[t]/eR)
			}else{
				#option 2, NOT expected
				
				## Guessing rAA[t] is not easy.. Blank just in case!
				logR[t]<- NaN
			}
		}
		
		# 6: Blank for BB while positive for AA & AB
		if( (e.tAA) & (e.tAB) & !(e.tBB)){
			#two options!
			# 1: SNP is "correctly between" tAA and tAB
			# 2: Heterozygous SNP is subjected to deletion or UPD of allele A making it unexectedly to be between tAB and tBB where it normally should not NOT BE.
				
			if(Theta[t]<= tAB[t]){
				#option 1! the expected
				#interpolate as SNP is expected to be between tAB and tBB
				eR<-rAA[t]+(Theta[t]-tAA[t])*(rAB[t]-rAA[t])/(tAB[t]-tAA[t])
				logR[t] <- log2(R[t]/eR)
			}else{
				#option 2, NOT expected
				## Guessing rBB[t] is not easy.. Blank just in case!
				logR[t]<- NaN
			}
		}
				
		# 7: positive for AA & BB & AB
		if((e.tAA) & (e.tAB) & (e.tBB)){
			# AA & BB & AB exists. DO Illumina style calculation
			if(Theta[t]<tAB[t]){
				#lower part of BAF scale (0-0.5 typically)
				#tAA<theta<tAB
				eR<-rAA[t]+(Theta[t]-tAA[t])*(rAB[t]-rAA[t])/(tAB[t]-tAA[t])
				logR[t] <- log2(R[t]/eR)
			}else{
				#upper part of BAF scale (0.5-1 typically)
				#tAB<theta<tBB
				eR<-rAB[t]+(Theta[t]-tAB[t])*(rBB[t]-rAB[t])/(tBB[t]-tAB[t])
				logR[t] <- log2(R[t]/eR)
			}
		}
	} #end for
	rm(tmp.tBB)
	rm(tmp.tAA)
	return(logR)
}
###


### Require the R package Limma from www.bioconductor.org
library(limma)
###

cat("Performing normalization\n")
## perform tQN normalization ###
for(r in 1:nbrAssays){
	#foreach master assay, calculate, write and plot for the 4 different platforms.
	sampleName<-assaysToRow$Assay[r]
	cat(paste("  - sample",sampleName," (",r,"/",nbrAssays,")\n"))
	output.file<-paste(sampleName,"_tQN.txt",sep="")
	output.file<-paste(output.path,output.file,sep="")
        if(!file.exists(output.file)){
        #### Read sample data ####
        #format (tab separated):
	#Name	Chr	Position	X	Y
          my.file<-paste(input.path,sampleName,"_extracted.txt",sep="")
          baf.data<-read.delim(my.file,header=TRUE,na.strings=c(NA,NaN))
          if(length(which(is.na(baf.data[,"X"])))/nrow(baf.data)>0.20){
            cat(paste("\t sample",sampleName,"has more than 20% missing X values => skipped\n"))
          }
          else{
            baf.data$Name<-as.character(baf.data$Name)
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

            print(head(AA))
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

            print(theta.tQN[1:10])
	#### Calculate tQN X and Y to fit theta and R
            Y.tQN<- R.tQN*tan(theta.tQN*pi/2) / (1+tan(theta.tQN*pi/2))
            X.tQN<-R.tQN-Y.tQN
	###
	
	### Add CNV probes back on if present ###
            if(cnv_present){
              theta.tQN<-c(theta.tQN,theta.cnv)
              R.tQN<-c(R.tQN,R.cnv)
              baf.data<-rbind(baf.data,cnv.data)
              X.tQN<-c(X.tQN,cnv.data$X)
              Y.tQN<-c(Y.tQN,cnv.data$Y)
              rm(cnv.data,R.cnv,theta.cnv)
              gc();
            }
	###
	
	## Calculate BAF from corrected Theta
            cat("      * finding matching probes\n")
            ii<-match(baf.data$Name,clusterFile.new$reporterId)
            dd<-cbind(baf.data,clusterFile.new[ii,-c(1)])
            if(length(which(is.na(ii)))>0){
		#NA exists!!
              dd<-dd[-which(is.na(ii)),]
              baf.data<-baf.data[-which(is.na(ii)),]
              theta.tQN<-theta.tQN[-which(is.na(ii))]
              R.tQN<-R.tQN[-which(is.na(ii))]
		#AA<-AA[-which(is.na(ii)),]
              X.tQN<-X.tQN[-which(is.na(ii))]
              Y.tQN<-Y.tQN[-which(is.na(ii))]
            }
            rm(ii)
            gc();

            cat("      * computing new baf\n")
	### BAF ###
            corrected.baf<-baf.function(theta.tQN,dd$AA_T_Mean,dd$AB_T_Mean,dd$BB_T_Mean)
	###
	
            cat("      * computing new LRR\n")

        ### Calculate Log R ###
            corrected.logR<-logR.function(R.tQN,theta.tQN,dd$AA_T_Mean,dd$AB_T_Mean,dd$BB_T_Mean,dd$AA_R_Mean,dd$AB_R_Mean,dd$BB_R_Mean)
	###
	
            rm(dd)
            gc();
	
	
	
	### Write out the corrected BAF and logRratio
            vv<-data.frame(Name=baf.data$Name,Chr=baf.data$Chr,Position=baf.data$Position,corrBAF=corrected.baf,corrLogRratio=corrected.logR,corrX=X.tQN,corrY=Y.tQN)
            output.file<-paste(sampleName,"_tQN.txt",sep="")
            output.file<-paste(output.path,output.file,sep="")
            write.table(vv, file =output.file, quote = FALSE,append=FALSE, sep = "\t", row.names = FALSE,col.names = c("Name","Chr","Position","tQN B Allele Frequency", "tQN Log R Ratio","tQN X","tQN Y"))
            rm(vv)
            gc();
	###
          }
        }
}
cat(paste("Analysis finished",date(),"\n\n"))

#end each master assay	
