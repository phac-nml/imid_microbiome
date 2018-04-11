

darkcols <- c(brewer.pal(8, "Accent"),rev(brewer.pal(8, "Dark2")[-8]), brewer.pal(8,"Paired"))
cbbPalette <- c( "#E69F00", "#56B4E9", "#8B2323","#66CD00","#9932CC", "#D55E00", "#00008B", "#CDC0B0")

if(!exists("colorchoice")){
	colorchoice=cbbPalette
}
#### input required depthThresh and presentThresh



### getting data from the object
###MRcounts(metaobj); pData(metaobj);fData(metaobj)

#Library sizes (depths of coverage) and normalization factors , for the samples
#libSize  normFactors



## modified on 2018/01 
if(!exists("patientid")){patientid="Patient"}
if(!exists("diseaseid")){diseaseid="Disease"}
if(!exists("sampleid")){sampleid="Sample_name"}
	

###############
###filter OTU presence or minimum depth 
###############
#present:Features with at least 'present' postive samples.
#depth:Samples with at least this much depth of coverage
##filterData(metaobj, present = 10, depth = 100)
#metaobj=filterData(metaobj, present = 10, depth = 1)

metaobj=filterData(metaobj, present = presentThresh, depth = depthThresh)
if(perct==""){
	resultDir=paste0("present",presentThresh,"_depth",depthThresh)
}else{
	resultDir=paste0("present",presentThresh,"_depth",depthThresh,"_perct",perct)
}

dir.create(resultDir)
pdf(paste0(resultDir,"/Hist_Sums_raw_counts_afterfilter_present",presentThresh,".pdf"))
hist(log2(colSums(MRcounts(metaobj))),xlab="log2(Sum of counts)", main="By Samples", breaks=50)
hist(log2(rowSums(MRcounts(metaobj))),xlab="log2(Sum of counts)", main="By OTUs", breaks=50)
dev.off()


###############
###Calculating normalization factors
###############
if(perct==""){
	perct = cumNormStat(metaobj) ##.53 , .5 w filter
	perct
}
metaobj = cumNorm(metaobj, p = perct)

###############
### way to export normalized count matrices and stats
###############
nmat = cbind(as.matrix(fData(metaobj)),MRcounts(metaobj, norm = TRUE, log = FALSE))
exportMat(nmat, file = paste0(resultDir,"/normedcnt_",perct,".tsv"))
rmat =  cbind(as.matrix(fData(metaobj)),MRcounts(metaobj, norm = FALSE, log = FALSE))
exportMat(rmat, file = paste0(resultDir,"/rawcnt_",perct,".tsv"))
exportStats(metaobj, p=perct,file = paste0(resultDir,"/normStats_",perct,".tsv"))
#head(read.csv(file = paste0(resultDir,"/normStats.tsv"), sep = "\t"))

pdf(paste0(resultDir,"/Hist_Sums_Normed_counts_afterfilter_present",presentThresh,".pdf"))
hist(log2(colSums(MRcounts(metaobj, norm = T))),xlab="log2(Sum of normed counts)", main="By Samples", breaks=50)
hist(log2(rowSums(MRcounts(metaobj, norm = T))),xlab="log2(Sum of normed counts)",, main="By OTUs", breaks=50)
dev.off()

### printing out the outlying samples in terms of total raw counts after filter
print(head(sort(log2(colSums(MRcounts(metaobj))))))
print(tail(sort(log2(colSums(MRcounts(metaobj))))))

pdf(paste0(resultDir,"/Boxplot_log2_counts_afterfilter_present",presentThresh,".pdf"), width=12)
boxplot( MRcounts(metaobj, norm = F, log = TRUE), main="Before Norm")
boxplot( MRcounts(metaobj, norm = TRUE, log = TRUE), main="After Norm")
dev.off()

###############
### data visualization
###############
table(pData(metaobj)[,diseaseid])
#CD HC MS RA UC 
#42 47 40 45 40 
#Crohn’s disease, ulcerative colitis, multiple sclerosis, rheumatoid arthritis) and healthy controls
#Each patient also had 2 stool samples (denoted as .1 or .2) 2 months apart. Biological replicates, add patientid as covariate
#r’s represent technical replicates. technical replicates, can use this to assess normalization!

#### data Heatmaps
require(gplots)
cl=factor(pData(metaobj)[,diseaseid])
#cols=brewer.pal(12, "Set3")
cols=colorchoice
clcol=cols[as.integer(cl)]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)


pdf(paste0(resultDir,"/Heatmap_sampleDistance_Bray-Curtis_logCnt_normB4Af_OTU.pdf"), width=9, height=9)
###Bray–Curtis dissimilarity
library(vegan)
sampBCTT=as.matrix(vegdist(t(MRcounts(metaobj, norm = T, log = T)), method="bray"))
heatmap.2(sampBCTT, trace="none", scale="none", ColSideColors=clcol, main="samp_Bray-Curtis,normT,logT", col=heatmapCols)
legend("topright",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
sampBCFT=as.matrix(vegdist(t(MRcounts(metaobj, norm = F, log = T)), method="bray"))
heatmap.2(sampBCFT, trace="none", scale="none", ColSideColors=clcol, main="samp_Bray-Curtis,normF,logT", col=heatmapCols)
dev.off()

pdf(paste0(resultDir,"/Heatmap_FeatureDistance_Bray-Curtis_logCnt_normB4Af_OTU.pdf"), width=9, height=9)
###Bray–Curtis dissimilarity
fcl=factor(fData(metaobj)$Class)
fcols=darkcols
fclcol=fcols[as.integer(fcl)]
library(vegan)
fBCTT=as.matrix(vegdist(MRcounts(metaobj, norm = T, log = T), method="bray"))
heatmap.2(fBCTT, trace="none", scale="none", ColSideColors=fclcol, main="samp_Bray-Curtis,normT,logT", col=heatmapCols)
legend("topright",  legend = levels(fcl), col = fcols[1:length(levels(fcl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
fBCFT=as.matrix(vegdist(MRcounts(metaobj, norm = F, log = T), method="bray"))
heatmap.2(fBCFT, trace="none", scale="none", ColSideColors=fclcol, main="samp_Bray-Curtis,normF,logT", col=heatmapCols)
dev.off()




### get the cor/BC between replicates (tech and bio)
jMeltAnnot<-function(sampCor){ ### note metaobj is used in the function

	require(reshape)
	sampCorTTm<-melt(sampCor)
	sampCorTTm<-sampCorTTm[sampCorTTm$X1!=sampCorTTm$X2,]
	sampCorTTm_1<-apply(sampCorTTm, 1,function(y)paste(sort(y[1:2]),collapse=","))
	sampCorTTm<-sampCorTTm[!duplicated(sampCorTTm_1),]

	addInfo<-t(apply( sampCorTTm,1,function(y){
		sameDis=pData(metaobj)[y[1],diseaseid]==pData(metaobj)[y[2],diseaseid]
		#samePat=pData(metaobj)[y[1],patientid]==pData(metaobj)[y[2],patientid]
		#techRep=gsub("r","",pData(metaobj)[y[1],"Sample_name"])==gsub("r","",pData(metaobj)[y[2],"Sample_name"])
		#bioRep=samePat& !techRep
		if(pData(metaobj)[y[1],patientid]==pData(metaobj)[y[2],patientid]){
			if(gsub("r","",pData(metaobj)[y[1],sampleid])==gsub("r","",pData(metaobj)[y[2],sampleid])){
				ret="techRep"}else{ret="sameIndiv"}
		}else{
			ret="diffPatients"
		}
		c(sameDis, ret)
	}))
	colnames(addInfo)=c("sameCategory", "patient")
	
	master=addInfo[,"patient"]
	master[addInfo[,"patient"]=="diffPatients"&addInfo[,"sameCategory"]==FALSE]="diffCategories"
	master[addInfo[,"patient"]=="diffPatients"&addInfo[,"sameCategory"]==TRUE]="sameCategory"
	master=factor(master, levels=c("techRep","sameIndiv","sameCategory","diffCategories"))
	data.frame(sampCorTTm,addInfo, master)
}

### plotting the index comparing techReps, before and after norm 
jPlotTechRepB4Aft<-function(mtableTT,mtableFT, distIndex="", resultDir="", alter="two.sided", feat="OTU"){
	if(resultDir!=""){
		resultDir=paste0(resultDir,"/")
	}
	
	if(sum(mtableTT$patient=="techRep")>1){
		mtableTTtech=mtableTT[mtableTT$patient=="techRep","value"]
		mtableFTtech=mtableFT[mtableFT$patient=="techRep","value"]
		
		print(quantile(mtableTTtech))
		print(quantile(mtableFTtech))
		print(paste("T.test",alter,"for",distIndex,"between techRep b4 and after norm",signif(t.test(mtableTTtech,mtableFTtech, alternative=alter)$p.value,4)))
		print(paste("wilcox.test",alter,"for",distIndex,"between techRep b4 and after norm",signif(wilcox.test(mtableTTtech,mtableFTtech, alternative=alter)$p.value,4)))
		print(length(mtableTTtech))
		pdf(paste0(resultDir,"TechAssessment_",distIndex,"_b4after_norm_",feat,"_.pdf"))
		 boxplot(cbind(normed=mtableTTtech,b4norm=mtableFTtech), main=paste(distIndex,"on techRep"))
		dev.off()
	}
	
	require(ggplot2)
	pdf(paste0(resultDir,"Density_",distIndex,"_samples_filtered_normed_",feat,"_2_.pdf"), width=9)
		print(ggplot(mtableTT,aes(x=value, fill=master)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex) + ggtitle("After Norm"))
		print(ggplot(mtableTT,aes(x=value, fill=patient)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex) + ggtitle("After Norm"))
		print(ggplot(mtableTT,aes(x=value, fill=sameCategory)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex)+ ggtitle("After Norm"))
		print(ggplot(mtableFT,aes(x=value, fill=master)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex)+ ggtitle("B4 Norm"))
		print(ggplot(mtableFT,aes(x=value, fill=patient)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex)+ ggtitle("B4 Norm"))
		print(ggplot(mtableFT,aes(x=value, fill=sameCategory)) + geom_density(alpha=0.25)+ xlim(0,1)+  xlab(distIndex)+ ggtitle("B4 Norm"))
	dev.off()
	
	require(ggplot2)
	pdf(paste0(resultDir,"Boxplot_",distIndex,"_samples_filtered_normed_",feat,"_2_.pdf"), width=9)
		c0=unique(mtableTT$master)#[c(4,3,2,1)] #c(4,1,2,3)
		ylim0=range(c(mtableTT$value,mtableFT$value))
		print(ylim0)
		boxplot(lapply(c0,function(y){mtableTT[mtableTT$master==y,"value"]}),ylim=ylim0, names=c0, main=paste(distIndex,"After Norm"))
		boxplot(lapply(c0,function(y){mtableFT[mtableFT$master==y,"value"]}),ylim=ylim0, names=c0, main=paste(distIndex,"Before Norm"))
	dev.off()
	
	require(ggpubr)
	#mtableTT$master<-factor(as.character(mtableTT$master),unique(mtableTT$master))
	#mtableFT$master<-factor(as.character(mtableFT$master),unique(mtableFT$master))
		#ggdotplot(mtableTT, "master", "value", fill = "master", color = "master",binwidth = 0.0025,  ylab = "Bray-Curtis distance", add = "boxplot")#median_iqr
	ggboxplot(mtableFT, "master", "value", #color = "master",  
		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
		ylab = "Bray-Curtis Dissimilarity", xlab="", ylim=c(0,0.9), main="Before Norm")%>%
	ggexport(filename = paste0(resultDir,"ggdotplot_",distIndex,"_samples_filtered_normed_",feat,"_b4norm.pdf"))

	ggboxplot(mtableTT, "master", "value", #color = "master",  
		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
		ylab = "Bray-Curtis Dissimilarity", xlab="", ylim=c(0,0.9), main="After Norm")%>%
	ggexport(filename = paste0(resultDir,"ggdotplot_",distIndex,"_samples_filtered_normed_",feat,"_Afnorm.pdf"))


}


### tech assessment: sample Bray-Curtis density plots
sampBCTTmOTU=jMeltAnnot(sampBCTT)
sampBCFTmOTU=jMeltAnnot(sampBCFT)
jPlotTechRepB4Aft(sampBCTTmOTU, sampBCFTmOTU, distIndex="BC", resultDir)
medianVals<-sapply(levels(sampBCTTmOTU$master), function(y){median(sampBCTTmOTU[sampBCTTmOTU$master==y,"value"])})
names(medianVals)=levels(sampBCTTmOTU$master)
IQRVals<-sapply(levels(sampBCTTmOTU$master), function(y){c(IQR=IQR(sampBCTTmOTU[sampBCTTmOTU$master==y,"value"]),quantile(sampBCTTmOTU[sampBCTTmOTU$master==y,"value"]))})
print("IQRvalues")
round(IQRVals,2)

### plotting the differences in distribution between diff categories and same
pdf(paste0(resultDir,"Boxplot_","BC","_samples_filtered_normed_OTU_diffSameCate.pdf"))
	a=sampBCTTmOTU[sampBCTTmOTU$patient=="diffPatients"&sampBCTTmOTU$sameCategory==T,"value"]
	b=sampBCTTmOTU[sampBCTTmOTU$patient=="diffPatients"&sampBCTTmOTU$sameCategory==F,"value"]
	wilcox.test(a,b,alter="less") #p-value < 2.2e-16
	boxplot(a,b,names=c("sameCategories","diffCategories"),ylab="Bray-Curtis Dissimilarity", main="After norm")
	c=sampBCFTmOTU[sampBCFTmOTU$patient=="diffPatients"&sampBCFTmOTU$sameCategory==T,"value"]
	d=sampBCFTmOTU[sampBCFTmOTU$patient=="diffPatients"&sampBCFTmOTU$sameCategory==F,"value"]
	wilcox.test(c,d,alter="less") #p-value < 2.2e-16
	boxplot(c,d,names=c("sameCategories","diffCategories"),ylab="Bray-Curtis Dissimilarity", main="Before norm")
dev.off()



#Aggregating counts
metaobjPhylF = aggTax(metaobj, lvl = "Phylum", out = "matrix",norm = T,log = F, aggfun = colSums)
head(metaobjPhylF[, 1:5])
metaobjClalF = aggTax(metaobj, lvl = "Class", out = "matrix",norm = T,log = F, aggfun = colSums)
head(metaobjClalF[, 1:5])
metaobjGen = aggTax(metaobj, lvl = "Genus", out = "matrix",norm =T,log = F, aggfun = colSums)
head(metaobjGen[, 1:5])
metaobjGenFF = aggTax(metaobj, lvl = "Genus", out = "matrix",norm =F,log = F, aggfun = colSums) #raw data


pdf(paste0(resultDir,"/Heatmap_logCnt_normed_OTUvsGenus.pdf"), width=9, height=9)
heatmap.2(MRcounts(metaobj, norm = T, log = T), trace="none", scale="none", ColSideColors=clcol, main="OTU,normT,logT", col=heatmapCols)
legend("topright",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
heatmap.2(MRcounts(metaobj, norm =F, log = T), trace="none", scale="none", ColSideColors=clcol, main="OTU,normF,logT", col=heatmapCols)
legend("topright",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
heatmap.2(log2(metaobjGen+1), trace="none", scale="none", ColSideColors=clcol, main="agg Genus, normedT,logged2+1", col=heatmapCols)
legend("topright",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
heatmap.2(log2(metaobjGenFF+1), trace="none", scale="none", ColSideColors=clcol, main="agg Genus, normedF,logged2+1", col=heatmapCols)
legend("topright",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=.7,box.col = NA)
dev.off()



sampBCTTmGenus=jMeltAnnot(as.matrix(vegdist(t(log2(metaobjGen+1)), method="bray")))
sampBCFTmGenus=jMeltAnnot(as.matrix(vegdist(t(log2(metaobjGenFF+1)), method="bray")))

jPlotTechRepB4Aft(sampBCTTmGenus,sampBCFTmGenus, distIndex="Bray-Curtis", resultDir, alter="two.sided", feat="Genus")



pdf(paste0(resultDir,"/Barplot_allSamples_relativeAbundance.pdf"), width=floor(sqrt(ncol(metaobjPhylF))))
bardat=metaobjPhylF; bardatInfo="Phylum, sum of counts"
bardat<-bardat[,order(pData(metaobj)[,diseaseid], pData(metaobj)[,patientid])]
par(mar=c(3,4,2,7))
 barplot(t(t(bardat)/colSums(bardat)*100), col=darkcols, ylab="Relative Abundance %", main=paste(bardatInfo,",ordered samples"),names.arg=rep("",ncol(bardat)))
 legend("topright", inset=c(-0.1,0),rownames(bardat),fill=darkcols[1:nrow(bardat)], xpd=T, cex=0.8, bty = "n")
bardat=metaobjClalF; bardatInfo="Class, sum of counts"
bardat<-bardat[,order(pData(metaobj)[,diseaseid], pData(metaobj)[,patientid])]
par(mar=c(3,4,2,7))
 barplot(t(t(bardat)/colSums(bardat)*100), col=darkcols, ylab="Relative Abundance %", main=paste(bardatInfo,",ordered samples"),names.arg=rep("",ncol(bardat)))
 legend("topright", inset=c(-0.12,0),rownames(bardat),fill=darkcols[1:nrow(bardat)], xpd=T, cex=0.8, bty = "n")
dev.off()



#### decisions!!! function to summarize to genus and function to summarize samples
require(robustbase)
metaobjClalFmed= aggTax(metaobj, lvl = "Class", out = "matrix",norm = T,log = F, aggfun = colMedians)
metaobjPhylFmed= aggTax(metaobj, lvl = "Phylum", out = "matrix",norm = T,log = F, aggfun = colMedians)
metaobjClalFmean= aggTax(metaobj, lvl = "Class", out = "matrix",norm = T,log = F, aggfun = colMeans)
metaobjPhylFmean= aggTax(metaobj, lvl = "Phylum", out = "matrix",norm = T,log = F, aggfun = colMeans)

jAggSamp<-function(y, xdata, group, fcn=median){
	apply(xdata[,group==y],1,fcn)
}

### current thought and from lit: sum (of OTUs belong to a phylum/class), then mean (of samples) is more robust
pdf(paste0(resultDir,"/Barplot_categories_relativeAbundance_differentWays_2.pdf"), width=9)
	udis=sort(unique(pData(metaobj)[,diseaseid]))
	phyCountsSumMean<-sapply(udis,jAggSamp, xdata=metaobjPhylF, group=pData(metaobj)[,diseaseid], fcn=mean)
	bardat=phyCountsSumMean; bardatInfo="Phylum, CountsSumMean"
	par(mar=c(3,3,2,7))
	 barplot(t(t(bardat)/colSums(bardat)*100), col=darkcols, main=bardatInfo)
	 legend("topright", inset=c(-0.25,0),rownames(bardat),fill=darkcols[1:nrow(bardat)], xpd=T, cex=0.8)

	claCountsSumMean<-sapply(udis,jAggSamp, xdata=metaobjClalF, group=pData(metaobj)[,diseaseid], fcn=mean)
	bardat=claCountsSumMean; bardatInfo="Class, CountsSumMean"
	par(mar=c(3,3,2,7))
	 barplot(t(t(bardat)/colSums(bardat)*100), col=darkcols, main=bardatInfo)
	 legend("topright", inset=c(-0.25,0),rownames(bardat),fill=darkcols[1:nrow(bardat)], xpd=T, cex=0.8)

dev.off()


####  MDS 
jMDS<-function(mydata,cl, titleAdd="", col2use=cbbPalette, save2pdf=F, plot3D=F, dirSaveIn=""){
	if(dirSaveIn!="") dirSaveIn=paste0(dirSaveIn,"/")
	mydata=t(mydata)
	# d <- dist(mydata) # euclidean distances between the rows 
	require(vegan)
	d <- vegdist(mydata, method="bray") # Bray-Curtis
	fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim
	#fit # view results

	# plot solution 
	x <- fit$points[,1]
	y <- fit$points[,2]
	if(save2pdf){
		pdf(paste0(dirSaveIn,"MDS_",titleAdd,"_Bray-Curtis.pdf"))
	}
	plot(x, y, xlab="MDS 1", ylab="Coordinate 2",  main=paste("Metric MDS", titleAdd), col=col2use[as.integer(cl)], pch=20)
	legend("topleft",levels(cl),col=col2use,pch=20)

	plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",  main=paste("Metric MDS", titleAdd),	type="n")
	text(x, y, labels = row.names(mydata), cex=.7, col=col2use[as.integer(cl)]) #col=clcol
	legend("topleft",levels(cl),col=col2use,pch=20)
	if(save2pdf){
		dev.off()
	}
	if(plot3D){
		require(rgl)
		z <- fit$points[,3]
		plot3d(x, y, z, col=col2use[as.integer(cl)], type="s", size =1, lwd=4)
		text3d(x,y,z+0.08, row.names(mydata), fontweight="bold", cex= 1, col=col2use[as.integer(cl)])

	}
}

cl = factor(classes)
save2pdf0=T
jMDS(mydata= MRcounts(metaobj, norm = F, log = TRUE), cl, titleAdd="normF_logT_OTU", col2use=colorchoice, save2pdf=save2pdf0, dirSaveIn=resultDir)
jMDS(mydata= MRcounts(metaobj, norm = TRUE, log = TRUE), cl,titleAdd= "normT_logT_OTU", col2use=colorchoice, save2pdf=save2pdf0, dirSaveIn=resultDir)
jMDS(mydata=log2(1+metaobjGen), cl, "normT_logT_Genus", col2use=colorchoice, save2pdf=save2pdf0, dirSaveIn=resultDir)

save.image(paste0(resultDir,"/initial_result.rdata"))

pdf(paste0(resultDir,"/Hist_OTU_Genus_count_bias.pdf"))
 hist(sort(table(fData(metaobj)$Genus)), breaks=30, xlab="Total numbers of OTUs per Genus", main="Total OTU counts per Genus")
dev.off()

print("done")

