########################################################################################################
# WGCNA
# 20211225
##################################
# first, do the quantile normalization
library(preprocessCore)
library('WGCNA')
library(cluster)
setwd('D:\\Projects\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcirptome_Xingxing')
met_2 <- read.csv('Metabolites_matrix_XX_exp1_overlapping_with_AS_expression.csv',row.names=1,head=T,stringsAsFactors=F)
#TPM_2 <- read.table('Exp_matrix_XX_exp1.txt',row.names=1,head=T,stringsAsFactors=F)
TPM_normalized=normalize.quantiles(as.matrix(TPM_2))
rownames(TPM_normalized) <- rownames(TPM_2)
# get rid of genes with TPM < 2 in all samples
TPM <- TPM_normalized[apply(TPM_normalized >= 2, 1, any),]
# further filter genes if there is no significant difference in any of the comparison
FC <-  read.csv('FC_Exp1_linking_to_XX_met_responsive_in_at_least_one_sample.csv',row.names=1,head=T,stringsAsFactors=F)
TPM <- TPM[rownames(TPM) %in% rownames(FC),]
colnames(TPM) <- colnames(TPM_2)
# save normalized TPM values, for genes which are expressed in at least one sample, and are differentially expressed in at least one cultivar in D2 vs D1
write.csv(TPM,'TPM_AS_Exp1_overlapping_with_XX_met_expressed_differential.csv',row.names=T,quote=F)

# save the median value across duplicates
library(matrixStats)
TPM_median <- c()
for(cul in c('VS16','WBC3','Exp1_F1')){
	subdat <- TPM[,grep(colnames(TPM),pattern=cul,fixed = TRUE)]
	for(days in c('D1','D2')){
		subtem <- subdat[,grep(colnames(subdat),pattern=days,fixed = TRUE)]
		TPM_median <- cbind(TPM_median,rowMedians(as.matrix(subtem)))
		colnames(TPM_median)[dim(TPM_median)[2]] = paste(cul,days,sep='_')
		}
	}
rownames(TPM_median) <- rownames(TPM)
TPM_median <- TPM_median[apply(TPM_median >= 2, 1, any),]
write.csv(TPM_median,'TPM_AS_Exp1_overlapping_with_XX_met_expressed_differential_median.csv',row.names=T,quote=F)

########################################################################################################
# WGCNA, https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
dat <- as.matrix(read.csv('TPM_AS_Exp1_overlapping_with_XX_met_expressed_differential_median.csv',row.names=1,head=T,stringsAsFactors=F))
dataExprVar <- dat
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# filter genes with small variation across samples
# m.mad <- apply(dataExpr,1,mad)
# dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
# soft threshold
## check for the outliers
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=70, by=2))
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('TPM_xingxing_Exp1_expressed_differential_median_soft.pdf')
par(mfrow = c(2,1))
cex1 = 0.9
#x-axis is Soft threshold (power)，y-axis is the scale free topology model fit，the higher the value on y-axis ,the closer to non-scale
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="blue")
abline(h=0.85,col="red")
	 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# R-square=0.9
dev.off()

# build network
enableWGCNAThreads()
softPower = 32 # change this value depending on what you get
adjacency = adjacency(dataExpr, power = softPower, type = "signed") #specify network type
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
library(flashClust)
geneTree = flashClust(as.dist(dissTOM), method="average")
pdf('TPM_AS_Exp1_overlapping_with_XX_met_Gene_Clustering_on_TOM-based_dissimilarity.pdf')
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(dataExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
pdf('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_Clustering_of_module_eigengenes.pdf')
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.1
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

# get the reordered genes in the tree, the order is the original order of genes
geneordered <-colnames(dataExpr)[geneTree$order]
# get the color of genes
colorlabel <- cbind(colnames(dataExpr),mergedColors)
rownames(colorlabel) <- colorlabel[,1]
color_order <- colorlabel[geneordered,2]
write.csv(color_order,'TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_colors.csv',row.names=T,quote=F)

#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
dev.off()

write.csv(MEs,'TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_MEs.csv',row.names=T,quote=F)

save(MEs, moduleLabels, moduleColors, geneTree, file= "TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_Network_allSamples_signed_nomerge_RLDfiltered.RData")


#################
# draw the z-score of expression for each module
# do the z-score transformation
library(rcompanion)
dat <- as.matrix(read.csv('TPM_AS_Exp1_overlapping_with_XX_met_expressed_differential_median.csv',row.names=1,head=T,stringsAsFactors=F))
zscore <- dat
for(i in 1:nrow(dat)){
	zscore[i,] <- blom(dat[i,],method='zscore')
	}
module <- read.csv('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_colors.csv',head=T,row.names=1)
zscore <- zscore[rownames(module),]
zscore <- cbind(zscore,module)
pdf('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_zscore_profile.pdf')
par(mfrow=c(3,3))
for(n in unique(module[,1])){
	subdat <- zscore[zscore[,ncol(zscore)]==n,1:(ncol(zscore)-1)]
	MAX <- max(subdat)
	MIN <- min(subdat)
	CN <- colnames(subdat)
	CNN <- length(CN) # number of samples
	# Plot the first gene
	plot(t(subdat[1,]),type="l",ylim=c(MIN,MAX),col="gray",xaxt="n",xlab="",ylab=n)
	axis(1,at=seq(1,CNN,length.out=ncol(subdat)),labels=CN,las=2)
	# Plot the rest of the genes
	for (i in seq(2,nrow(subdat))){
	  lines(t(subdat[i,]),col="gray")
	}
	# Plot the median
	MED <- sapply(subdat,median)
	lines(MED,col="red")
	}
dev.off()

# draw the correlation among modules
MEs_col = MEs
#colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# correlation among modules according to gene expression
# marDendro/marHeatmap
pdf('Correlation_among_AS_exp1_overlapping_with_XX_met_WGCNA_modules_median_MEDissThres_0.1.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
#########################
# associate the expression modules with traits
trait <- "Metabolites_matrix_XX_exp1_overlapping_with_AS_expression_median.txt"
traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
colnames(traitData) <- rownames(dataExpr)						
traitData <- t(traitData)

### association between modules and phenotype data 
nSamples <- nrow(dataExpr)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('Correlation_between_XX_exp1_metabolites_and_WGCNA_modules_median_MEDissThres_0.1.pdf')
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
write.csv(modTraitCor,'Correlation_between_XX_exp1_metabolites_and_WGCNA_modules_median_MEDissThres_0.1.csv',row.names=T,quote=F)
Per95 <- quantile(modTraitCor[!is.na(modTraitCor)],0.95) #0.8004502

####################################################################################
########## 																	########
##########                      WGCNA of metabolites
########## 																	########
####################################################################################
dat <- as.matrix(read.table('Metabolites_matrix_XX_exp1_overlapping_with_AS_expression_median.txt',row.names=1,head=T,stringsAsFactors=F,sep='\t'))
dat <- dat[apply(dat > 0, 1, any),]
dataExprVar <- dat
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# filter genes with small variation across samples
#m.mad <- apply(dataExpr,1,mad)
#dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
# soft threshold
## 
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=70, by=2))
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('Metabolites_XX_Exp1_overlapping_with_AS_expression_soft.pdf')
par(mfrow = c(2,1))
cex1 = 0.9
# 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="blue")
abline(h=0.85,col="red")
	 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#
dev.off()

# build network
enableWGCNAThreads()
softPower = 16
adjacency = adjacency(dataExpr, power = softPower, type = "signed") #specify network type
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method="average")
library(flashClust)
geneTree = flashClust(as.dist(dissTOM), method="average")
pdf('Metabolites_XX_Exp1_overlapping_with_AS_expression_Gene_Clustering_on_TOM-based_dissimilarity.pdf')
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(dataExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEs[is.na(MEs)] <- 0
MEDiss= 1-cor(MEs)
MEDiss[is.na(MEDiss)] <- 0
METree= flashClust(as.dist(MEDiss), method= "average")
pdf('Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_Clustering_of_module_eigengenes.pdf')
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.1
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

# get the reordered genes in the tree, the order is the original order of genes
geneordered <-colnames(dataExpr)[geneTree$order]
# get the color of genes
colorlabel <- cbind(colnames(dataExpr),mergedColors)
rownames(colorlabel) <- colorlabel[,1]
color_order <- colorlabel[geneordered,2]
write.csv(color_order,'Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_colors.csv',row.names=T,quote=F)

#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
dev.off()

write.csv(MEs,'Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_MEs.csv',row.names=T,quote=F)
 
#save(MEs, moduleLabels, moduleColors, geneTree, file= "Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_Network_allSamples_signed_nomerge_RLDfiltered.RData")


#################
# draw the z-score of expression for each module
# do the z-score transformation
library(rcompanion)
dat <- as.matrix(read.table('Metabolites_matrix_XX_exp1_overlapping_with_AS_expression_median.txt',row.names=1,head=T,stringsAsFactors=F,sep='\t'))
dat <- dat[apply(dat > 0, 1, any),]
zscore <- dat
for(i in 1:nrow(dat)){
	zscore[i,] <- blom(dat[i,],method='zscore')
	}
module <- read.csv('Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_colors.csv',head=T,row.names=1)
zscore <- zscore[rownames(module),]
zscore <- cbind(zscore,module)
pdf('Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_zscore_profile.pdf')
par(mfrow=c(3,3))
for(n in unique(module[,1])){
	subdat <- zscore[zscore[,ncol(zscore)]==n,1:(ncol(zscore)-1)]
	MAX <- max(subdat)
	MIN <- min(subdat)
	CN <- colnames(subdat)
	CNN <- length(CN) # number of samples
	# Plot the first gene
	plot(t(subdat[1,]),type="l",ylim=c(MIN,MAX),col="gray",xaxt="n",xlab="",ylab=n)
	axis(1,at=seq(1,CNN,length.out=ncol(subdat)),labels=CN,las=2)
	# Plot the rest of the genes
	for (i in seq(2,nrow(subdat))){
	  lines(t(subdat[i,]),col="gray")
	}
	# Plot the median
	MED <- sapply(subdat,median)
	lines(MED,col="red")
	}
dev.off()


##################################################################
# draw the correlation among metabolites and modules expression modules
MEs <- read.csv('Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_MEs.csv',head=T,row.names=1)
ME_expression <- read.csv('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_MEs.csv',head=T,row.names=1)
colnames(ME_expression) <- paste('TPM_',colnames(ME_expression),sep='')
ME_both <- cbind(ME_expression,MEs)
MEs_col = ME_both
#MEs_col = orderMEs(MEs_col)
pdf('Correlation_among_expression_and_metabolites_modules_XX_exp1.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

ME_expression <- read.csv('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_MEs.csv',head=T,row.names=1)
MEs <- read.csv('Metabolites_XX_Exp1_overlapping_with_AS_expression_WGCNA_MEDissThres_0.1_MEs.csv',head=T,row.names=1)
colnames(MEs) <- paste('Metabolites_',colnames(MEs),sep='')
ME_both <- cbind(ME_expression,MEs)
MEs_col = ME_both
#MEs_col = orderMEs(MEs_col)
pdf('Correlation_among_expression_and_metabolites_modules_XX_exp1_02.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()




# draw the correlation among modules
MEs_col = MEs
#colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# marDendro/marHeatmap 
pdf('Correlation_among_AS_exp1_overlapping_with_XX_met_WGCNA_modules_median_MEDissThres_0.1.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

