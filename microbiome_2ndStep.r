### using out of bag stats

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Go in the posData or allData folder, double click on "metaobj_b4filter.rdata" to load
##### run in bulk w Jess_microbiome_prelim_Batch.r
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 	

require(reshape)
require(metagenomeSeq)
require(ggpubr)
require(Rcpp)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(robustbase)

### criteria to remove samples w low depth 
depthThresh=1
### criteria to remove OTUs w low presence among samples
#presentThresh<-floor(length(pData(metaobj)$Disease)/length(unique(pData(metaobj)$Disease)))
presentThresh<-40
perct=0.75


jesslabels=c("CD", "HC", "MS", "RA", "UC", "Diseased")
jessCol=c("purple1", "goldenrod1", "springgreen1", "orangered1", "dodgerblue1", "darkgrey")
colorchoice=jessCol
source("../microbiome_prelim_Batch.r")
	
	
	
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### double click into the folder on "initial_result.rdata" to load
##### run in bulk w microbiome_machineLearning_Batch.r
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 	
	
require(randomForest)
require(caret)
require(metagenomeSeq)
require(ROCR)
require(e1071)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

	
#perc4Val=0.3  ## the percentage of left out data for validation
ntree=500
saveInFold<-"OOB"

############### all pair-wise
positives=c('CD','MS','RA','UC','HC')
comparisons2<-NULL
for (iii in 1:(length(positives)-1)){
	for(jjj in 2:length(positives)){
		if(jjj>iii)
			comparisons2<-c(comparisons2,paste(positives[iii],positives[jjj], sep="2"))
	}
}
comparisons2<-c(comparisons2,"Diseased")

dataiss=c("OTU", "Genus")
bAccMat<-matrix(0,nrow=length(comparisons2),ncol=2)
rownames(bAccMat)=comparisons2
colnames(bAccMat)=dataiss

BAccuracyList0<-list()
meanprecRec0List<-list()
bAcc0Matmean<-bAccMat
bAcc0MatmeanRep<-bAccMat
###aucV0 and aucV
aucV0mean<-bAccMat
aucV0meanRep<-bAccMat
varImpMergeMeanL<-list(OTU=NULL, Genus=NULL)



for(compIs in comparisons2){
	if(compIs!='Diseased'){
		positive=strsplit(compIs,"2")[[1]][1]
	}else{
		positive=compIs
	}
	for(datai in dataiss){

		print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		load("initial_result.rdata")	
		
		datis=datai;samp2exlude="bio"; repis=1
		
		source("../../microbiome_machineLearning_Batch.r")
			print(datai)
		
		bAcc0Matmean[compIs,datai]=meanBAccuracy0
		
		aucV0mean[compIs,datai]=mean(aucV0)
		varImpMergeMeanL[[datai]]<-cbind(varImpMergeMeanL[[datai]],  varImpMergeMean)
		colnames(varImpMergeMeanL[[datai]])[ncol(varImpMergeMeanL[[datai]])]=compIs
	
		#meanprecRec0List[[compIs]][[datai]]=meanprecRec0
		meanprecRec0List[[datai]]=rbind(meanprecRec0List[[datai]],meanprecRec0)
		BAccuracyList0[[compIs]][[datai]]=BAccuracy0

		bAcc0MatmeanRep[compIs,datai]=repBAmean0
		aucV0meanRep[compIs,datai]=mean(aucV0Rep)
		
	}
}
varImpMergeMeanL
bAcc0Matmean
aucV0mean
meanprecRec0List
meanprecRec0List2<-lapply(meanprecRec0List, function(y){rownames(y)=rownames(bAcc0Matmean);y})

save.image("pair-wise_OOB.rdata")


##############################
### save the variable importance in an Excel file
##############################
require(openxlsx)
hs <- createStyle(textDecoration = "BOLD", fontSize=11, fontName="Arial") #, fgFill = "#4F80BD", fontColour = "#FFFFFF"
### also requires Rtools
Sys.setenv("R_ZIPCMD" = "C:\\Rtools\\bin\\zip.exe")  ### Rtools was downloaded and installed on the machine, temp set of envir
excelName<-paste0(saveInFold,paste("/varImp",samp2exlude,"ntree",ntree,"_annot.xlsx", sep="_"))
#write.xlsx(DEtables, file = excelName)
print(paste("Saving varImp results in .xlsx file:", excelName))
#varImpMergeMedList<-c(varImpMergeMedianRMVF, varImpMergeMedianRMVT)
#names(varImpMergeMedList)=paste(names(varImpMergeMedList),c(rep("RMVF",2),rep("RMVT",2)), sep="_")
varImpMergeMeanL2<-lapply(varImpMergeMeanL,function(y) data.frame(ID=rownames(y),y))
## adding annot
varImpMergeMeanL2$OTU= cbind(as.matrix(fData(metaobj)),varImpMergeMeanL2$OTU)
genus2OTU=sapply(as.character(varImpMergeMeanL2$Genus[,"ID"]), function(y){tid=fData(metaobj)$Genus==y;paste(fData(metaobj)[tid,"Group"],collapse=",")})
varImpMergeMeanL2$Genus=cbind(genus2OTU,varImpMergeMeanL2$Genus)
write.xlsx(varImpMergeMeanL2, file = excelName, headerStyle = hs)


##############################
### save the performance in an Excel file
##############################
excelName<-paste0(saveInFold,paste("/OOB_Performance",samp2exlude,"ntree",ntree,".xlsx", sep="_"))
#write.xlsx(DEtables, file = excelName)
print(paste("Saving OOB_Performance results in .xlsx file:", excelName))
perfList<-c(list(BA=bAcc0Matmean, AUC=aucV0mean, BArep=bAcc0MatmeanRep,AUCrep=aucV0meanRep),meanprecRec0List2)
perfList<-lapply(perfList,function(y) data.frame(ID=rownames(y),y))
write.xlsx(perfList, file = excelName, headerStyle = hs)





#######################################################################
##### multi-class
#######################################################################

ntree=500
saveInFold<-"OOB"

positive='cate';compIs=positive
dataiss=c("OTU", "Genus")
#overallaccM=NULL
overallaccM=list()
meanBAacrossclass=overallaccM
varImpMergeMeanLM=list()
for(datai in dataiss){

	print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print(datai)
	load("initial_result.rdata")	
	
	datis=datai;samp2exlude="bio"; repis=1
	
	source("../../microbiome_machineLearning_Batch.r")
		
	mean( rfResults[[1]]$confusionm$byClass[,"Balanced Accuracy"])
	meanBAacrossclass[[datai]]= sapply(rfResults, function(y){mean(y$confusionm$byClass[,"Balanced Accuracy"])})
	
	#overallaccM[[datai]]=c(overallaccM[[datai]],overallacc)
	overallaccM[[datai]]=c(overallaccM[[datai]],overallacc)
	

	varImpMergeMeanLM[[datai]]<-varImpMergeMean

}
save(positive,compIs, samp2exlude,repis,overallaccM, varImpMergeMeanLM, file="MultiGroup_OOB.rdata")
mean( meanBAacrossclass$OTU)# 0.6869427  #all dat 0.6889682
mean( meanBAacrossclass$Genus) #0.6697899 # all dat  0.6666268






