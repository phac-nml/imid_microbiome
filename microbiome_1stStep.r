
require(metagenomeSeq)
setwd("C:\\Users\\jchen\\_NML\\microbiome_codeSummary")
#############
### load data
#############

### load count and taxa, 1st column name was added as OTU from Jess' file
#dataTab<-read.table("otus_tax_julie_11sept2017.txt", header=T, sep="\t", as.is=T)
dataTab<-read.table("otumat_taxmat_forjulie_21dec2017.txt", header=T, sep="\t", as.is=T)
dataTab<-dataTab[,c(6,1:5,7:ncol(dataTab))]

###count
cntdat<-list(counts=dataTab[,7:ncol(dataTab)], taxa=dataTab[,"Group"])
rownames(cntdat$counts)=cntdat$taxa
###taxa
taxadat<-dataTab[,1:6]
rownames(taxadat)=taxadat[,1]

### load metadata
phenodat<-read.table("sample_data_3.txt", header=T, sep="\t", as.is=T)
rownames(phenodat)=paste0("X", phenodat$SampleID)


### X1204.2.1 in count file is X1204.2r in metadata file
colnames(cntdat$counts)[colnames(cntdat$counts)!=rownames(phenodat)]
# [1] "X1204.2.1" "X1235.1.1" "X1257.2.1" "X2011.2.1" "X2019.2.1" "X2093.2.1" "X2196.1.1" "X3074.1.1" "X3295.2.1" "X5011.2.1" "X5013.2.1" "X8124.1.1" "X8170.2.1"
rownames(phenodat)[colnames(cntdat$counts)!=rownames(phenodat)]
# [1] "X1204.2r" "X1235.1r" "X1257.2r" "X2011.2r" "X2019.2r" "X2093.2r" "X2196.1r" "X3074.1r" "X3295.2r" "X5011.2r" "X5013.2r" "X8124.1r" "X8170.2r"

rownames(phenodat)=gsub("r$",".R",rownames(phenodat))
colnames(cntdat$counts)=rownames(phenodat)


###############
#### setup the MRexperiment object
###############
phenotypeData = AnnotatedDataFrame(phenodat)
phenotypeData

OTUdata = AnnotatedDataFrame(taxadat)
OTUdata

metaobj = newMRexperiment(cntdat$counts,phenoData=phenotypeData,featureData=OTUdata)

pdf("Hist_Sums_raw_counts_new2018.pdf")
hist(log2(colSums(MRcounts(metaobj))),xlab="log2(Sum of counts)", main="By Samples", breaks=50)
hist(log2(rowSums(MRcounts(metaobj))),xlab="log2(Sum of counts)",, main="By OTUs", breaks=50)
dev.off()

dir.create("allData")
save.image("allData/metaobj_b4filter.rdata")


################################################################################
#### only gram positive
########################################################################
print(nrow(cntdat$counts))#7943

gramPos<-!taxadat$Phylum%in%c("Verrucomicrobia","Proteobacteria","Bacteroidetes","Bacteria_unclassified","Candidatus_Saccharibacteria","Chloroflexi", "Fusobacteria","Planctomycetes","Spirochaetes","Synergistetes","Lentisphaerae", "Deinococcus-Thermus")
#Deinococcus-Thermus, have thick cell walls that give them gram-positive stains, but they include a second membrane and so are closer in structure to those of gram-negative bacteria.
#The Firmicutes are a phylum of bacteria, most of which have Gram-positive cell wall structure. A few, however, such as Megasphaera, Pectinatus, Selenomonas and Zymophilus, have a porous pseudo-outer membrane that causes them to stain Gram-negative
###gram positive and negative table https://en.wikipedia.org/wiki/Tenericutes
### pos: "Firmicutes"          "Actinobacteria"    "Tenericutes"   

cntdat$counts=cntdat$counts[gramPos,]
cntdat$taxa=cntdat$taxa[gramPos]

###taxa
taxadat<-taxadat[gramPos,]
print(nrow(cntdat$counts)) #7943 -> 7118

###############
#### setup the MRexperiment object
###############
phenotypeData = AnnotatedDataFrame(phenodat)
phenotypeData

OTUdata = AnnotatedDataFrame(taxadat)
OTUdata

metaobj = newMRexperiment(cntdat$counts,phenoData=phenotypeData,featureData=OTUdata)

pdf("Hist_Sums_raw_counts_new2018_PosOnly.pdf")
hist(log2(colSums(MRcounts(metaobj))),xlab="log2(Sum of counts)", main="By Samples", breaks=50)
hist(log2(rowSums(MRcounts(metaobj))),xlab="log2(Sum of counts)",, main="By OTUs", breaks=50)
dev.off()

hist(rowSums(MRcounts(metaobj)!=0),xlab="Sum of samples w nonzero counts",, main="By OTUs", breaks=50)

dir.create("posData")
save.image("posData/metaobj_b4filter.rdata")
