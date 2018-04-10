### 2018 0208 added balanced sampling 
library(randomForest)
library(caret)
library(metagenomeSeq)
library(ROCR)

####
### inputs
####

#datis="Genus"
#datis="OTU"
#samp2exlude="bio"; repis=1

##### multi-class comparison
#positive='cate';compIs=positive
##### two-class comparison
#positive='CD';compIs="CD2HC"
##### two-class comparison
#positive='Diseased';compIs=positive


#perc4Val=0.3  ## the percentage of left out data for validation
#ntree=500
 
####
###scripts
####

perc4Val=NULL ## no left out
#saveInFold<-"OOB"
if( !dir.exists(saveInFold)){
	dir.create(saveInFold)
}
### bias of using OTU

datclass=pData(metaobj)$Disease
if(datis=="OTU"){
	dat4model=MRcounts(metaobj, norm = T, log = TRUE)#; datis="OTU"
	print(paste("All",nrow(dat4model)))

}
if(datis=="Genus"){
	dat4model=log2(1+metaobjGen)#; datis="Genus"
	print(paste("All",nrow(dat4model)))

}

dat4model00<-dat4model


### samples to exclude
if(samp2exlude=="tech"){
	###excluding .r technical replicates for prediction, keep in mind there are still bio replicates too
	techrid=grep("r",pData(metaobj)$Sample_name)
	dat4model<-data.frame(cate=datclass[-techrid],t(dat4model[,-techrid]))
}
if(samp2exlude=="bio"){
	### exclude bio &tech replicates too.
	if(repis==1){
		repid=grep("\\.2|\\.2r|\\.1r",pData(metaobj)$Sample_name)
	}else{
		repid=grep("\\.1|\\.2r|\\.1r",pData(metaobj)$Sample_name)
	}
	dat4model<-data.frame(cate=datclass[-repid],t(dat4model[,-repid]))
}
dat4model0=dat4model




if(compIs=="cate"){
	### sample 4 validation
	print(as.matrix(table(dat4model$cate)))
	## get id for each category to sample
	catids<-lapply(unique(dat4model$cate),function(y){which(dat4model$cate==y)})
	names(catids)=unique(dat4model$cate)
}else if(length(grep("2",compIs))>0){
	#positive='CD'
	cate2learn=strsplit(compIs, "2")[[1]]
		### using only CD and normal
		dat4model<-dat4model0[dat4model0$cate%in%cate2learn,]
		dat4model$cate<-droplevels(dat4model$cate)
		
		### sample 4 validation
	print(as.matrix(table(dat4model$cate)))
	## get id for each category to sample
	catids<-lapply(unique(dat4model$cate),function(y){which(dat4model$cate==y)})
	names(catids)=unique(dat4model$cate)
}else if(compIs=="Diseased"){
	dat4model<-dat4model0
	dat4modelcate0=dat4model$cate
	dat4model$cate=as.character(dat4model$cate)
	### using only Diseased and normal
	dat4model$cate[dat4model$cate!="HC"]="Diseased"
	dat4model$cate=factor(dat4model$cate)

	### sample 4 validation
	print(as.matrix(table(dat4model$cate)))
	## get id for each category to sample
	catids<-lapply(unique(dat4model$cate),function(y){which(dat4model$cate==y)})
	names(catids)=unique(dat4model$cate)
} 


########################	 
### setting up training and validation data, no long doing so wihout left out
########################
###for each item in list: 	 dat4train, dat4vald
if( !is.null(perc4Val)){
	seedNums=1:10
	datTV<-lapply(seedNums, function(seedNum,perc4Val){

		set.seed(seedNum)
		valdids<-sapply(catids,function(y){
			ycnt2samp<-floor(length(y)*perc4Val)
			sample(y,ycnt2samp)
		})
		
		dat4vald<-dat4model[sort(unlist(valdids)),]
		dat4train<-dat4model[-sort(unlist(valdids)),]	
		
		list(dat4train=dat4train, dat4vald=dat4vald)
	},perc4Val=perc4Val)
 
 }

 
########################	 
### model building with the list of training and validation datasets
######################## 
dat4train=dat4model
seedNums=1:10

### for unbalanced data; 2 fold or more counts, sample balanced
sampCnt<-as.vector(table(dat4train[,1]))
if(length(sampCnt)==2 & abs(log2(sampCnt[1]/sampCnt[2]))>=1){
	sampCnt<-rep(min(sampCnt),2)
}

rfResults<-lapply(seedNums, function(seedNum){
	set.seed(seedNum)
	#Fit Random Forest Model
	rf1 = randomForest(cate ~ .,  
					   ntree = ntree,
					   sampsize=sampCnt,
					   data = dat4train)
	print(rf1) 

	# Variable Importance
	varImpPlot(rf1,   sort = T, n.var=10, main="Top 10 - Variable Importance")
	var.imp = data.frame(importance(rf1,type=2))
	var.imp$Variables = row.names(var.imp)  
	#print(var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),])
	
	### this is from OOB result
	print(cmv<-confusionMatrix(data= rf1$predicted,  
					reference=dat4train$cate,
					positive=positive))
						
						
	list(model=rf1, var.imp=var.imp, confusionm=cmv)
})
print("randomForestDone")


varImpMerge=do.call(cbind,lapply(rfResults,function(y) y$var.imp[,1]))
rownames(varImpMerge)=rfResults[[1]]$var.imp[,2]
varImpMergeMedian=rowMedians(varImpMerge) ## not too far off from median
names(varImpMergeMedian)=rownames(varImpMerge)
varImpMergeMean=rowMeans(varImpMerge) ## not too far off from median
names(varImpMergeMean)=rownames(varImpMerge)

write.table(varImpMerge,
	file=paste0(saveInFold,paste("/varImpMerge",datis,samp2exlude,compIs,".txt", sep="_")),
	col.names=T, row.names=T, quote=F,sep="\t")




##################
### NEW!!! with ComplexHeatmap: heatmap and bar plot
##################
plottop=10
varImpMerge2<-data.frame(mean0=apply(varImpMerge,1,mean),se0=apply(varImpMerge,1,sd)/sqrt(ncol(varImpMerge)))
varImpMerge2<-varImpMerge2[order(varImpMerge2$mean0, decreasing=T),] ## order of varImpMerge2 diff from varImpMerge
	   
col2use=c(colorchoice[1:6],colorchoice[2])
names(col2use)[1:7]=c("CD","HC","MS","RA","UC", "Diseased", "Healthy")
require(gplots)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
pdf(paste0(saveInFold,paste("/ComplexHeatmap_dataFromTopVarImp",datis,samp2exlude,compIs,"_horizLegend.pdf", sep="_")), width=11,height=7)

	hdat2plot=t(data.matrix(dat4model[,head(rownames(varImpMerge2), plottop)]))
	colnames(hdat2plot)=rep("   ", ncol(hdat2plot))
	#ha2 = rowAnnotation( p2 = row_anno_barplot(colMeans(t(varImpMerge[head(rownames(varImpMerge2), plottop),])), axis = TRUE),  width = unit(5, "cm"))
	ha2 = rowAnnotation( b2 = row_anno_boxplot(varImpMerge[head(rownames(varImpMerge2), plottop),], axis = TRUE, outline=F),  width = unit(3, "cm"))
	ha1 = HeatmapAnnotation(df = data.frame(cohort=as.character(dat4model[,1])), show_legend = TRUE, col = list(cohort = col2use), annotation_legend_param = list(title="Cohort",legend_direction="horizontal", nrow=1))

	colCH=brewer.pal(9,"PuBuGn")
	h2plot<-Heatmap(hdat2plot,cluster_rows =F,col=colCH, top_annotation = ha1, top_annotation_height = unit(1, "cm"), heatmap_legend_param = list(title = "Abundance", color_bar = "continuous", legend_direction="horizontal")) + ha2
	print(h2plot)
dev.off()

pdf(paste0(saveInFold,paste("/ComplexHeatmap_dataFromTopVarImp",datis,samp2exlude,compIs,".pdf", sep="_")), width=9,height=7)

	hdat2plot=t(data.matrix(dat4model[,head(rownames(varImpMerge2), plottop)]))
	colnames(hdat2plot)=rep("   ", ncol(hdat2plot))
	#ha2 = rowAnnotation( p2 = row_anno_barplot(colMeans(t(varImpMerge[head(rownames(varImpMerge2), plottop),])), axis = TRUE),  width = unit(5, "cm"))
	ha2 = rowAnnotation( b2 = row_anno_boxplot(varImpMerge[head(rownames(varImpMerge2), plottop),], axis = TRUE, outline=F),  width = unit(5, "cm"))
	ha1 = HeatmapAnnotation(df = data.frame(cohort=as.character(dat4model[,1])), show_legend = TRUE, col = list(cohort = col2use), annotation_legend_param = list(title="Cohort", annotation_legend_side = "bottom"))

	colCH=brewer.pal(9,"PuBuGn")
	h2plot<-Heatmap(hdat2plot,cluster_rows =F,col=colCH, top_annotation = ha1, top_annotation_height = unit(1, "cm"), heatmap_legend_param = list(title = "Abundance", color_bar = "continuous")) + ha2
	print(h2plot)
dev.off()

			   
### for 2 class comparisons, save ROC curves
if(length(levels(dat4train[,1]))==2){
	pdf(paste0(saveInFold,paste("/ROCR",datis,samp2exlude,compIs,"ntree",ntree,".pdf", sep="_")), width=5, height=5)

		#jBalAcc<-function(mat0){mat0<-mat0[,1:nrow(mat0)];mean(diag(mat0)/rowSums(mat0))}
		#BAccuracy0<-sapply(rfResults,function(y) jBalAcc(y$model$confusion))
		
		BAccuracy0<-sapply(rfResults,function(y) y$confusionm$byClass["Balanced Accuracy"])
		
		medBAccuracy0=median(BAccuracy0) ## balanced
		meanBAccuracy0=mean(BAccuracy0) ## balanced
		print(BAccuracy0)
		print(medBAccuracy0)
		print(meanBAccuracy0)
		precRec0<-sapply(rfResults,function(y) y$confusionm$byClass[c("Precision","Recall")])
		meanprecRec0<-rowSums(precRec0)/length(seedNums)
		
		aucV0<-rep(0,length(seedNums))
		colors <- 1:10
		for (i in 1:length(seedNums)) {
			rf.pred0<-rfResults[[i]][[1]]$votes[,positive]
			trueCate0<-dat4train$cate
			y.te0=as.numeric(trueCate0==positive)
			predROC0 <- prediction(rf.pred0,y.te0)
			aucV0[i]=signif(performance(predROC0, measure = "auc")@y.values[[1]],2)
				
			plot(performance(predROC0, measure = "tpr", x.measure = "fpr"), add=(i!=1),col=colors[i],lwd=2, lty=i, main=paste0("OOB ",paste(levels(trueCate0),collapse=" vs "),": ", datis, " meanBalAcc: ",signif(medBAccuracy0,2)))
			abline(a=0, b= 1, col="darkgray",lty=2)	
			}
		text(0.7,0.2,labels=paste("meanAUC=",mean(aucV0)))
		
		##odd that they can be slightly different lengths...
		##performance(predROC0, measure = "tpr", x.measure = "fpr")@y.values[[1]])

	
	
		if(FALSE){ #one seed only
			####
			#### from cross validation, training data that was OOB
			####
			### my balanced accuracy function: NOTE, the rows have to be the TRUE counts, while columns are prediction, otherwise it's computed the wrong way
			jBalAcc<-function(mat0){mat0<-mat0[,1:nrow(mat0)];mean(diag(mat0)/rowSums(mat0))}
			
			BAccuracy0<-jBalAcc(rf1$confusion)
			print(BAccuracy0)
			
			
			library(ROCR)
			rf.pred0<-rf1$votes[,positive]
			trueCate0<-dat4train$cate
			y.te0=as.numeric(trueCate0==positive)
			predROC0 <- prediction(rf.pred0,y.te0)
			aucV0=signif(performance(predROC0, measure = "auc")@y.values[[1]],2)
					
			plot(performance(predROC0, measure = "tpr", x.measure = "fpr"),col="darkgray",lwd=2, lty=i, main=paste0("OOB ",paste(levels(trueCate0),collapse=" vs "),": ", datis, " BalAccuracy: ",signif(BAccuracy0,2)))
			abline(a=0, b= 1, col="darkgray",lty=2)	

			text(0.8,0.2,labels=paste("AUC.OOB=",mean(aucV0)))
		}
	
	dev.off()
	
}

			   
				   
				   
				   
### for multi-class comparisons				   
if(length(levels(dat4train$cate))>2){
	print("Overall accuracies")
	overallacc<-sapply(rfResults,function(y) y$confusionm$overall["Accuracy"])
	print(overallacc)
	print(paste("Mean overall accuracy",mean(overallacc)))
	### OTU for multiclass: 0.4741; Genus for multiclass: 0.4370
	
	lapply(rfResults,function(y) y$confusionm$byClass)

	rowMedians(sapply(rfResults,function(y) y$confusionm$byClass[,"Sensitivity"]))
	rowMedians(sapply(rfResults,function(y) y$confusionm$byClass[,"Specificity"]))
	rowMedians(sapply(rfResults,function(y) y$confusionm$byClass[,"Precision"]))
	rowMedians(sapply(rfResults,function(y) y$confusionm$byClass[,"Recall"]))
	#print(rowMedians(sapply(rfResults,function(y) y$confusionm$byClass[,"Balanced Accuracy"])))
}

				   


####################
## predict on the other set 2 months after
####################
#######################################
### prediction on the other set (.1 or .2)
## training on .1 and predicting on .2 yielded  0.9365 balanced accuracy (OTU, CD vs HC)
#######################################
print("prediction on data 2 month later")
if(samp2exlude=="bio"){

	### exclude bio &tech replicates too.
	if(repis==2){
		repid=grep("\\.2|\\.2r|\\.1r",pData(metaobj)$Sample_name)
	}else{
		repid=grep("\\.1|\\.2r|\\.1r",pData(metaobj)$Sample_name)
	}
	
	
	if(compIs=="cate"){
		dat4modelrep<-data.frame(cate=datclass[-repid],t(dat4model00[,-repid]))

	}else if(length(grep("2",compIs))>0){
		dat4modelrep<-data.frame(cate=datclass[-repid],t(dat4model00[,-repid]))
		dat4modelrep<-dat4modelrep[which(dat4modelrep$cate%in%names(catids)),]
		dat4modelrep$cate<-droplevels(dat4modelrep$cate)
	
	}else if(compIs=="Diseased"){
		dat4modelrep<-data.frame(cate=datclass[-repid],t(dat4model00[,-repid]))
		dat4modelrep$cate=as.character(dat4modelrep$cate)
		
		### using only Diseased and normal
		dat4modelrep$cate[dat4modelrep$cate!="HC"]="Diseased"
		dat4modelrep$cate=factor(dat4modelrep$cate)
	} 

}

### sample 4 validation
print(as.matrix(table(dat4model$cate)))
catids<-lapply(unique(dat4model$cate),function(y){which(dat4model$cate==y)})


jPredNew<-function(mmodel,newDat, positive){

	# Predicting response variable
	newDat$predicted.response <- predict(mmodel ,newDat[,-1])
	predProb<-predict(mmodel ,newDat[,-1], type = "prob")
	# Create Confusion Matrix
	print(  
	cmv<-confusionMatrix(data=newDat$predicted.response,  
					reference=newDat$cate,
					positive=positive))
	list(dat4vald=newDat,predProb=predProb, confusionm=cmv)

}
repPred<-jPredNew(rfResults[[1]]$model,newDat=dat4modelrep, positive=positive)
rfResultsRep<-lapply(rfResults, function(y){
	jPredNew(y$model,newDat=dat4modelrep, positive=positive)
})

#repBAmean0<-mean(sapply(rfResultsRep,function(y) y$confusionm$byClass["Balanced Accuracy"])) 


if(length(levels(rfResultsRep[[1]]$dat4vald$cate))==2){
	pdf(paste0(saveInFold,paste("/ROCR",datis,samp2exlude,compIs,"ntree",ntree,"_rep.pdf", sep="_")), width=5, height=5)

			#medAccuracy=median(sapply(rfResults,function(y) y$confusionm$overall["Accuracy"]))  
			BAccuracy0Rep<-sapply(rfResultsRep,function(y) y$confusionm$byClass["Balanced Accuracy"])
			repBAmean0=mean(BAccuracy0Rep) ## balanced
			print(repBAmean0)
							
			library("ROCR")
			n <- length(rfResultsRep) #  n models
			colors <- 1:n #
			aucV0Rep<-rep(0,length(n))
			for (i in 1:n) {
				rf.pred=rfResultsRep[[i]]$predProb[,positive]
				y.te=as.numeric(rfResultsRep[[i]]$dat4vald$cate==positive)
				predROC <- prediction(rf.pred,y.te)
				aucV0Rep[i]=signif(performance(predROC, measure = "auc")@y.values[[1]],2)
				
				 plot(performance(predROC, measure = "tpr", x.measure = "fpr"), add=(i!=1),col=colors[i],lwd=2, lty=i, main=paste0(paste(levels(rfResults[[1]]$dat4vald$cate),collapse=" vs "),": ", datis, " medBalAccuracy: ",signif(repBAmean0,2)))
				abline(a=0, b= 1, col="darkgray",lty=2)
			
			}
			text(0.8,0.2,labels=paste("meanAUC=",mean(aucV0Rep)))
			
			
	
	dev.off()
}





print("done")


