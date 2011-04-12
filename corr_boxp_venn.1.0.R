#!/usr/bin/env Rscript

#Uses read counts as indata
#Plots and calculates R^2 values for data.frame object. Saves linear regression plots to disc with 
#column names as filename.
#Plots boxplots unsing Spearman correlation
#Draws Venn diagrams within replicates and across replicates. 
#Run with "Rscript linreg_script.R [infile.txt] [SampleName1] [SampleName2] 
#[Cond1_nr_of_col] [Cond2_nr_of_col] [Expressed_read_count_mean] [out_directory]" from bash
#©2011 Henrik Stranneheim

#1st argument infile
#2nd argument sample name 1
#3nd argument sample name 2
#4th argument condition 1 number of columns
#5th argument condition 2 number of columns
#6th argument expressed read counts median
#7th argument out directory

args <- commandArgs(TRUE) #Collects arguments

infile <-args[1] #Collect infile

sampleName1 <- args[2] #Collects sample name 1

sampleName2 <- args[3] #Collects sample name 2

nrcol_T <- as.numeric(args[4]) #Collects condition 1 nr of columns. NOTE: 1st column in file is always GeneName

nrcol_N <- as.numeric(args[5]) #Collects condition 2 nr of columns. NOTE: 1st column in file is always GeneName

expRCmedian <- as.numeric(args[6]) #The cut-off for counting a gene as expressed in a condition

od <- args[7] #Collect location of out directory

library(gplots) #For Venn diagrams

#Read infile
RC <-read.table(paste(infile), header=T)

#Select columns condition 1
RC_T<-RC[1:nrcol_T+1] #Adjusting for first column is gene name
#Add column names
colNameRC_T <-names(RC_T)

#Select column condition 2
RC_N<-RC[(2+nrcol_T):(1+nrcol_T+nrcol_N)]
#Add column names
colNameRC_N <-names(RC_N)

#set working directory
setwd( paste(od, sep="") )

#Printing to summmary file
#cat("Arguments:", "\n", "Working directory: ", "\t", getwd(), "/", "\n", sep="",
#file="corr_boxp_venn.1.0_out.txt", append=FALSE)

cat("Arguments:", "\n","Infile:", "\t", infile, "\n", "Out directory:", "\t", od, "\n",
"SampleName1:", "\t", sampleName1, "\n", "Nr of Columns:", "\t", nrcol_T, "\n",
"SampleName2:", "\t", sampleName2, "\n", "Nr of Columns:", "\t", nrcol_N, "\n",
"Expressed Gene if Read Count Median >", "\t", expRCmedian,"\n", sep="",
file="corr_boxp_venn.1.0_out.txt", append=FALSE)

#Printing to summmary file
cat(sampleName1, "columns:", "\t", colNameRC_T, "\n",file="corr_boxp_venn.1.0_out.txt", append=TRUE)
cat(sampleName2, "columns:", "\t", colNameRC_N, "\n", "\n",file="corr_boxp_venn.1.0_out.txt", append=TRUE)


#VennDiagram A431

system("mkdir -p VennDiagram")

setwd("./VennDiagram")

RC_Tmedian <- apply(RC_T,  1, median)
RC_Nmedian <- apply(RC_N,  1, median)

prT<-(RC_Tmedian>=expRCmedian)
prN<-(RC_Nmedian>=expRCmedian)

input_T_N.df <-data.frame(prT, prN)
names(input_T_N.df)[1] <- sampleName1
names(input_T_N.df)[2] <- sampleName2

pdf( paste( paste( sampleName1, "VS", sampleName2, "median", "VennD", sep="_"), "pdf", sep=".") )
venn(input_T_N.df)
dev.off()

#Writes if expressed or not
write.table(input_T_N.df,paste( paste( sampleName1, sampleName2, "boolean", "expressed", sep="_"),"txt", sep="."), sep="\t" )


#Determines unique expression for conditions and corresponding GeneName
t <- 0
n<-0
nt<-0
#unique_prT <-0
#unique_prN <-0
shared<-0
unRCT<-0
unRCTmedian<-0
unRCN<-0
unRCNmedian<-0
for(i in 1: length(prT)  )  {

	if ( prT [i] != prN[i] ) {
		
		if (prT [i] == "TRUE") {
		t <-t+1
		unRCT[t]<-as.character(RC[i,1])
		unRCTmedian[t]<-RC_Tmedian[i]
		
		}
		if (prN [i] == "TRUE") {
		n <-n+1
		unRCN[n]<-as.character(RC[i,1])
		unRCNmedian[n]<-RC_Nmedian[i]

		}
	}
	else  {
		nt<-nt+1
		#shared[nt] <- as.character(RC[i,1])
	}
}

##Plotting the unique entries median read counts

#Extracts the maximum value (among other things) 
sum1<-summary(unRCTmedian)
sum2<-summary(unRCNmedian)

#Printing to summmary file
setwd("..")

cat(sampleName1, "\n", "Nr of unique genes:", "\t", length(unRCT), "\n", 
"Mean Read Counts", "\t", sum1[4], "\n", "Median Read Counts", "\t", sum1[3], "\n", "\n", sep="", file="corr_boxp_venn.1.0_out.txt", append=TRUE)
cat(sampleName2, "\n", "Nr of unique genes:", "\t", length(unRCN), "\n", 
"Mean Read Counts", "\t", sum2[4], "\n", "Median Read Counts", "\t", sum2[3], "\n", "\n", sep="", file="corr_boxp_venn.1.0_out.txt", append=TRUE)

setwd("./VennDiagram")


if (length(unRCT) < length(unRCN) ) { #Determines the x-scale
	
	if ( sum1[6]<sum2[6]) { #Determines the y-scale
		
		pdf( paste( paste( sampleName1, "VS", sampleName2, "Unique", "RCmedian", sep="_"), "pdf", sep=".") )
		plot(c(1, length(unRCN) ), c(1,sum2[6]), col="white", xlab="Nr of Genes", ylab="Median Read Counts" )
		points(unRCTmedian, col="black", pch=21, lty=2)
		points(unRCNmedian, col="red", pch=22)
		legend(x="topright", c(paste( sampleName1),paste( sampleName2), 
		paste( paste(sampleName1, "All Genes RC Median", sep="-"), sum1[3], sep="=" ),
		paste( paste(sampleName2, "All Genes RC Median", sep="-"), sum2[3], sep="=" ) ),
		cex=0.8, col=c("black","red"), pch=21:24 )
		dev.off()

	}	else {
		
			pdf( paste( paste( sampleName1, "VS", sampleName2, "Unique", "RCmedian", sep="_"), "pdf", sep=".") )
			plot(c(1, length(unRCN) ), c(1,sum1[6]), col="white", xlab="Nr of Genes", ylab="Median Read Counts" )
			points(unRCTmedian, col="black", pch=21, lty=2)
			points(unRCNmedian, col="red", pch=22)
			legend(x="topright", c(paste( sampleName1),paste( sampleName2), 
			paste( paste(sampleName1, "All Genes RC Median", sep="-"), sum1[3], sep="=" ),
			paste( paste(sampleName2, "All Genes RC Median", sep="-"), sum2[3], sep="=" ) ),
			cex=0.8, col=c("black","red"), pch=21:24 )
			dev.off()
		}
}	else {

		if ( sum1[6]<sum2[6]) {
		pdf( paste( paste( sampleName1, "VS", sampleName2, "Unique", "RCmedian", sep="_"), "pdf", sep=".") )
		plot(c(1, length(unRCT) ), c(1,sum2[6]), col="white", xlab="Nr of Genes", ylab="Median Read Counts" )
		points(unRCTmedian, col="black", pch=21, lty=2)
		points(unRCNmedian, col="red", pch=22)
		legend(x="topright", c(paste( sampleName1),paste( sampleName2), 
		paste( paste(sampleName1, "All Genes RC Median", sep="-"), sum1[3], sep="=" ),
		paste( paste(sampleName2, "All Genes RC Median", sep="-"), sum2[3], sep="=" ) ),
		cex=0.8, col=c("black","red"), pch=c(21:24) )
		dev.off()
		
		} else {
			
			pdf( paste( paste( sampleName1, "VS", sampleName2, "Unique", "RCmedian", sep="_"), "pdf", sep=".") )
			plot(c(1, length(unRCT) ), c(1,sum1[6]), col="white", xlab="Nr of Genes", ylab="Median Read Counts" )
			points(unRCTmedian, col="black", pch=21, lty=2)
			points(unRCNmedian, col="red", pch=22)
			legend(x="topright", c(paste( sampleName1),paste( sampleName2), 
			paste( paste(sampleName1, "All Genes RC Median", sep="-"), sum1[3], sep="=" ),
			paste( paste(sampleName2, "All Genes RC Median", sep="-"), sum2[3], sep="=" ) ),
			cex=0.8, col=c("black","red"), pch=21:24 )
			dev.off()
			
			}
	}
	


#Creates data frame
unRCT.df <-data.frame(unRCT)
unRCN.df <-data.frame(unRCN)
names(unRCT.df)[1] <- sampleName1
names(unRCN.df)[1] <- sampleName2

#Writes corresponding geneName
write.table(unRCT.df, paste( paste( sampleName1, "Uniquely", "expressed", sep="_"),"txt", sep="."), sep="\t" )

write.table(unRCN.df, paste( paste( sampleName2, "Uniquely", "expressed", sep="_"),"txt", sep="."), sep="\t" )


#Venn diagram within replicates; condition 1

if ( (nrcol_T) < 5) {

vennNames <- c("A","B","C","D")
tot_T.df <- data.frame(c(1:nrow(RC_T) ) )

	for(i in 1: (nrcol_T)  )  {

		tot_T.df[i] <- (RC_T[,i]>=expRCmedian)
		names(tot_T.df)[i] <- vennNames[i]
	}
	pdf( paste( paste( sampleName1, "VennD", sep="_"), "pdf", sep=".") )
	venn(tot_T.df)
	dev.off()
} else {
	
	setwd("..")
	#Printing to summmary file
	cat(sampleName1,":", "\t", "VennDiagram cannot be drawn since the number of replicates are more than 4", "\n"
	,file="corr_boxp_venn.1.0_out.txt", append=TRUE)
	setwd("./VennDiagram")
	}
	
#Venn diagram within replicates; condition 2

if ( (nrcol_N) < 5) {

#vennNames <- c("A","B","C","D")
tot_N.df <- data.frame(c(1:nrow(RC_N) ) )

	for(i in 1: (nrcol_N)  )  {

		tot_N.df[i] <- (RC_N[,i]>=expRCmedian)
		names(tot_N.df)[i] <- vennNames[i]
	}
	pdf( paste( paste( sampleName2, "VennD", sep="_"), "pdf", sep=".") )
	venn(tot_N.df)
	dev.off()
} else {
	
	setwd("..")
	#Printing to summmary file
	cat(sampleName2, ":", "\t", "VennDiagram cannot be drawn since the number of replicates are more than 4", "\n"
	,file="corr_boxp_venn.1.0_out.txt", append=TRUE)
	setwd("./VennDiagram")
	}

######################
##### Regression ##### 
######################

setwd("..")

system("mkdir -p Regression")

setwd("./Regression")

#Condition 1

#Create variables
r2_T <- 0
z<-0

#For all columns except self
for(i in 1:( length(RC_T) -1 ) ) {
 for(k in (1+i):( (length(RC_T)) ) ) {
	z<-z+1 #Counter for R^2 value
	pdf(paste(colNameRC_T[i], colNameRC_T[k], "pdf", sep=".") )
	plot(RC_T[,i] ~ RC_T[,k], xlab=paste(colNameRC_T[i]), ylab=paste(colNameRC_T[k]))
	reg <-lm(RC_T[, i] ~RC_T[, k])
	abline(reg)
	legend("topleft", bty="n", legend=paste("R^2 = ",format(summary(reg)$adj.r.squared, digits=4)))
	dev.off()
	r2_T[z] <- as.numeric( format(summary(reg)$adj.r.squared, digits=4) )

	}
}

#Condition 1

#Create variables
r2_N<-0
z<-0

#For all columns except self
for(i in 1:( length(RC_N) -1 ) ) {
 for(k in (1+i):( (length(RC_N)) ) ) {
	z<-z+1 #Counter for R^2 value
	pdf(paste(colNameRC_N[i], colNameRC_N[k], "pdf", sep=".") )
	plot(RC_N[,i] ~ RC_N[,k], xlab=paste(colNameRC_N[i]), ylab=paste(colNameRC_N[k]))
	reg <-lm(RC_N[, i] ~RC_N[, k])
	abline(reg)
	legend("topleft", bty="n", legend=paste("R^2 = ",format(summary(reg)$adj.r.squared, digits=4)))
	dev.off()
	r2_N[z] <- as.numeric( format(summary(reg)$adj.r.squared, digits=4) )

	}
}


######################
#Spearman correlation#
######################
RCNoName <- RC[ , -1 ] #Removing gene names
x <- diag( nrow=length(RCNoName) ) #Creates a diagonal matrix

for(i in 1:length(RCNoName) ) { #Columns 

	for(j in 1:length(RCNoName)) { #Rows in x, columns in RCNoName
	
		x[i,j] <-cor(RCNoName[,i],RCNoName[,j], method="spearman", use="pairwise.complete.obs")
	}
}

#Selecting spearman correlation for condition 1
k <-0
sp_r2_T <-0
for(i in 1:( length(RC_T) -1) ) { #Row

	for(j in (1+i):length(RC_T) ) { #Column
		
		k <- k+1
		sp_r2_T[k] <- as.numeric( x[j,i] )
	}
}
#Selecting spearman correlation for condition 2
k <-0
sp_r2_N <-0
for(i in length(RCNoName):(length(RC_T) +2) ) { #Row

	for(j in (i-1):( length(RC_T)+1 ) ) { #Column
		
		k <- k+1
		sp_r2_N[k] <- as.numeric( x[j,i] )
	}
}


######################
#Creates vectors of equal length (required for data frame)
######################

if (length(r2_T) < length(r2_N) ) {

		r2_T[length(r2_T)+1:( length(r2_N)-length(r2_T) )] <-NA #Lin. reg
		sp_r2_T[length(sp_r2_T)+1:( length(sp_r2_N)-length(sp_r2_T) )] <-NA #Sperman. reg
		
}
if (length(r2_T) > length(r2_N)){
		
		r2_N[length(r2_N)+1:( length(r2_T)-length(r2_N) )] <-NA
		sp_r2_N[length(sp_r2_N)+1:( length(sp_r2_T)-length(sp_r2_N) )] <-NA
}

#Linear regression
r2<-as.data.frame(r2_T, row.names = c(1:length(r2_T)) )
r2[, 2]<-as.data.frame(r2_N)
names(r2)[1] <- sampleName1
names(r2)[2] <- sampleName2

#Spearman
sp_r2<-as.data.frame(sp_r2_T, row.names = c(1:length(sp_r2_T)) )
sp_r2[, 2]<-as.data.frame(sp_r2_N)

names(sp_r2)[1] <- sampleName1
names(sp_r2)[2] <- sampleName2

######################
#Boxplots#############
######################

#Linear regression
write.table(r2,paste( paste( sampleName1, sampleName2,"R_square_values", sep="_"),"txt", sep="."), sep="\t" )

pdf(paste( paste( sampleName1, sampleName2,"Boxplot", sep="_"),"pdf", sep=".") )
boxplot(r2_T,r2_N, 
names=c(paste( sampleName1), paste( sampleName2 ) ) )
title("R^2 values per condition")
dev.off()

#Spearman
write.table(x,paste( paste( sampleName1, sampleName2,"R_square_values_all_spearman", sep="_"),"txt", sep="."), sep="\t" )
write.table(sp_r2,paste( paste( sampleName1, sampleName2,"R_square_values_spearman", sep="_"),"txt", sep="."), sep="\t" )

pdf(paste( paste( sampleName1, sampleName2,"Boxplot_spearman", sep="_"),"pdf", sep=".") )
boxplot(sp_r2_T,sp_r2_N, 
names=c(paste( sampleName1), paste( sampleName2 ) ) )
title("R^2 values per condition (Spearman)")
dev.off()


