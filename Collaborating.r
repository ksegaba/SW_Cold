# Note for discussion
	1. For Xingxing exp data, exp1 had 2518 metabolites, while exp2 had only 733 metabolites. Then should we do the correlation analysis between metabolome and transcriptome separately for Exp1 and Exp2?
	2. Currently I did the analysis using RNA-seq data and Metabolites data for Exp1, both are from Xingxing. If we want to use RNA-seq data from Anne-Sophie, need to double check whether the sample ID between Anne-Sophie and Xingxing are the same. 
	
# cold responsive genes in Arabidopsis
CBF1	AT4G25490
CBF2	AT4G25470
CBF3	AT4G25480
DDF1	AT1G12610
MED2	AT1G11760	Pavir.3KG090754.1.p
ICE1	AT3G26744	Pavir.5KG747200.1.p, Pavir.5NG613400.1.p
ICE2	AT1G12860	Pavir.5KG747200.1.p, Pavir.5NG613400.1.p
CRPK1	AT1G16670	Pavir.7NG398000.1.p, Pavir.7KG349700.2.p
COR78	AT5G52310
FZF	AT2G24500	Pavir.4NG283800.1.p, Pavir.4KG386200.1.p
HOS1	AT2G39810	Pavir.9NG111000.1.p, Pavir.9KG518000.2.p
TCF1	AT3G55580	Pavir.5KG625700.6.p, Pavir.5NG587400.3.p
GLH2	AT4G04920	Pavir.9NG255973.1.p, Pavir.9KG137700.1.p
SVK		AT4G07395	
CPL1	AT4G21670	Pavir.7NG336900.1.p, Pavir.7KG258800.1.p, Pavir.1KG407900.4.p, Pavir.1NG367600.1.p
ZAT12	AT5G59820
B1L	AT1G18740	Pavir.7NG418000.1.p, Pavir.7KG363000.1.p, Pavir.5KG032100.2.p, Pavir.3NG176825.1.p, Pavir.9KG547500.1.p, Pavir.8KG043751.1.p, Pavir.8NG045100.1.p, Pavir.3KG039300.1.p, Pavir.3NG025600.1.p

LOS1	AT1G56070	Pavir.5KG523607.1.p, Pavir.5KG527914.1.p, Pavir.7NG380600.4.p, Pavir.5NG490100.1.p
LOS2	AT2G36530	Pavir.4NG033600.2.p, Pavir.4NG032400.1.p, Pavir.4KG019700.1.p
PER3	AT1G05260	Pavir.4NG293100.1.p, Pavir.4KG403700.2.p
U2BL	AT1G06960	Pavir.9KG471900.1.p, Pavir.9NG710000.2.p
GWD	AT1G10760	Pavir.4KG123300.1.p, Pavir.4NG147400.1.p
GNOM	AT1G13980	Pavir.9NG170800.1.p, Pavir.9KG063200.1.p
COR47	AT1G20440
LTI29	AT1G20450
LTI30	AT3G50970
ATL80	AT1G20823	Pavir.1KG452200.1.p, Pavir.1NG429100.1.p
ACC1	AT1G36160	Pavir.7KG040600.1.p, Pavir.7NG025300.2.p, Pavir.9NG401600.2.p, Pavir.9KG357200.2.p
COR413IM2	AT1G29390	Pavir.3NG181379.1.p, Pavir.3KG309200.1.p
COR413IM1	AT1G29395	Pavir.3NG181379.1.p, Pavir.3KG309200.1.p
AtRZ-1b	AT1G60650	Pavir.9KG049628.1.p, Pavir.9NG158500.1.p
ATRZ-1A	AT3G26420	Pavir.9KG169161.1.p, Pavir.9NG199400.4.p
AtRZ-1c	AT5G04280	Pavir.9KG049628.1.p, Pavir.9NG158500.1.p
PTP1	AT1G71860	Pavir.3KG051100.1.p, Pavir.3NG032300.3.p
AHB1	AT2G16060	Pavir.9KG538300.2.p, Pavir.9NG842400.1.p
CSP3	AT2G17870	
ATCSP4	AT2G21060	Pavir.1NG049200.3.p, Pavir.1KG014586.1.p
CAMTA3	AT2G22300	Pavir.9KG645000.2.p, Pavir.9NG818000.1.p, Pavir.9KG253013.1.p, Pavir.9KG253200.2.p, Pavir.9NG356600.2.p
AT2G23680	AT2G23680	

1NG367600	darkolivegreen
7NG398000	darkolivegreen
4KG386200	darkolivegreen
8NG045100	darkolivegreen
7KG349700	antiquewhite4
3KG039300	antiquewhite4
3NG025600	antiquewhite4
9KG253200	antiquewhite4
9KG518000	lightcoral
7NG380600	lightcoral
4NG147400	lightcoral
9NG199400	lightcoral
1KG014586	lightcoral
5KG625700	purple
3NG032300	purple
5NG587400	salmon4
4NG033600	salmon4
9KG357200	salmon4
5KG527914	honeydew
5NG490100	honeydew
4KG123300	darkviolet
9NG401600	blue4
9KG169161	coral1
9NG818000	darkseagreen4
9KG253013	darkseagreen2

#######################################################################################
###############															###############
###############			GO enrichment for WGCNA modules 20220216		###############
##############															###############
#######################################################################################
setwd('D:\\Projects\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcirptome_Xingxing\\Exp1\\GO_enrichment')
dat <- read.csv('TPM_AS_Exp1_overlapping_with_XX_met_WGCNA_MEDissThres_0.1_colors.csv',head=T)
dat[,1] <- gsub("Pavir.","",dat[,1])
dat[,1] <- gsub(".v5.1","",dat[,1])

for(colo in unique(dat$x)){
	subdat <- dat[dat$x == colo,]
	write.table(subdat[,1],paste(colo,'_gene_list.txt',sep=''),quote=F,sep='\t',row.names=F,col.names=F)
	}


############### GO enrichment
Enrichment <- function(k,n,C,G){
	return((k / C) / (n/ G))
	}
library('GOstats')
library('GSEABase')
library(pvclust)
library(gplots)
GO <- read.table('D:\\Projects\\Transformer_switchgrass\\05_GO_enrichment\\Switchgrass_GO_annotation.txt',head=F,sep='\t')
null <- rownames(read.table('D:\\Projects\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcirptome_Xingxing\\Exp1\\Exp1_read_counts_unique.txt',head=T,row.names=1))
null <- gsub("Pavir.","",null)
null <- gsub(".v5.1","",null)
null_go <- unique(merge(GO,t(t(null)),by.x='V1',by.y='V1'))
null_num <- length(null)

GO_pipeline <- function(Genes_file,short_name_to_save){
	pos <- read.table(Genes_file,head=F,sep='\t')
	pos[,1] <- gsub("Pavir.","",pos[,1])
	pos[,1] <- gsub(".v5.1","",pos[,1])
	pos_go <- unique(merge(GO,pos,by.x='V1',by.y='V1'))
	pos_num <- nrow(pos)
	res <- c()
	go <- unique(GO[,2])
	for(i in 1:length(go)){
		Pos <- pos_go[pos_go[,2]==go[i],]
		Null <- null_go[null_go[,2]==go[i],]
		a = nrow(Pos) # genes in cluster and with GO term
		b = pos_num - a # genes in cluster but with no GO term
		cc = nrow(Null) - a # genes with GO but not in cluster
		d = null_num - a - b - cc# genes with no GO and also not in cluster
		out <- c(go[i],a,b,cc,d)
		if(a+cc != 0) res <- rbind(res,out)
		}
	colnames(res) <- c('GO','Genes_in_module_with_GO','Genes_in_module_with_no_GO','Genes_not_in_module_with_GO','Genes_not_in_module_with_no_GO')
	write.table(res,paste('GO_',short_name_to_save,'.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	enrichment <- c()
	for(i in 1:nrow(res)){
		numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
		p <- fisher.test(numbers, alternative = "two.sided")[[1]][1]
		a = as.numeric(res[i,2])
		b = as.numeric(res[i,3])
		cc = as.numeric(res[i,4])
		d = as.numeric(res[i,5])
		if(Enrichment(a, a+b, a+cc, a+b+cc+d) >= 1) direction = '+' else direction = '-'
		enrichment <- rbind(enrichment, c(res[i,],direction,p))
		}
	write.table(enrichment,paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- read.table(paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
	dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
	dat <- dat[order(dat[,8]),]  
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- cbind(dat,'BP'='')
	dat <- cbind(dat,'CC'='')
	dat <- cbind(dat,'MF'='')
	for(i in 1:nrow(dat)){
		tryCatch( {
			if(!is.null(getGOTerm(dat[i,1])$BP[1])) dat[i,9] <- getGOTerm(dat[i,1])$BP[1]
			if(!is.null(getGOTerm(dat[i,1])$CC[1])) dat[i,10] <- getGOTerm(dat[i,1])$CC[1]
			if(!is.null(getGOTerm(dat[i,1])$MF[1])) dat[i,11] <- getGOTerm(dat[i,1])$MF[1]
			},
			error = function(e) {print(paste("no GO for ",dat[i,1]));NaN},
			finally = {})
		}
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	
	subdat <- dat[dat[,8] < 0.05 & dat$BP != '',]
	subdat <- cbind(subdat,'logP'=0)
	for(i in 1:nrow(subdat)){
		if(subdat[i,6] == '-')subdat[i,ncol(subdat)] = log10(subdat[i,8])
		if(subdat[i,6] == '+')subdat[i,ncol(subdat)] = -log10(subdat[i,8])
		}
	subdat[subdat$logP < -10,12] <- -10
	subdat[subdat$logP > 10,12] <- 10
	df <- as.matrix(cbind(subdat$logP,c(10,-10)))
	rownames(df) <- subdat$BP
	df <- df[order(df[,1],decreasing=T),]
	pdf(paste("GO_", short_name_to_save,'.fisher.qvalue_GO_term_BP.pdf',sep=''))
	heatmap.2(df,
		trace="none", # remove trace lines
		col = colorRampPalette(c("blue","white","red"))(21), # color palette
		dendrogram="none", # remove dendrogram
		cexRow=0.8,# row label font size
		labCol=F, # remove column labels
		Rowv=F, # turn off row clustering
		Colv=F, # turn off column clustering
		lmat=rbind(4:3,2:1),
		lwid=c(0.8,4), # dimensions of display array cell widths
		lhei=c(0.8,4)) # dimensions of display cell heights
	dev.off()
	write.table(subdat[,c(9,8,ncol(subdat))],paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term_BP.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')

	}


colo = 'darkseagreen2'
GO_pipeline(paste(colo,'_gene_list.txt',sep=''),colo)


1NG367600	darkolivegreen
7NG398000	darkolivegreen
4KG386200	darkolivegreen
8NG045100	darkolivegreen
7KG349700	antiquewhite4
3KG039300	antiquewhite4
3NG025600	antiquewhite4
9KG253200	antiquewhite4
9KG518000	lightcoral
7NG380600	lightcoral
4NG147400	lightcoral
9NG199400	lightcoral
1KG014586	lightcoral
5KG625700	purple
3NG032300	purple
5NG587400	salmon4
4NG033600	salmon4
9KG357200	salmon4
5KG527914	honeydew
5NG490100	honeydew
4KG123300	darkviolet
9NG401600	blue4
9KG169161	coral1
9NG818000	darkseagreen4
9KG253013	darkseagreen2









#######################################################################################
###############															###############
###############			fold change 20210216							###############
##############				betweem days								###############
#######################################################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# BiocManager::install(c("edgeR", "limma"))
library("edgeR")
library("limma")
library('stringr')
setwd('D:\\Projects\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcirptome_Xingxing')
# dat_AP13 <- read.table('AP13_read_counts_unique.txt',head=T,sep='\t',row.names=1)
# dat_exp2 <- read.table('Exp2_read_counts_unique.txt',head=T,sep='\t',row.names=1)
counts_all = read.table("Exp1_read_counts_unique.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
counts <- counts_all[,grep(colnames(counts_all),pattern='AS',fixed = TRUE)]
met <- read.table("Metabolites_matrix_XX_exp1.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
colnames(met) <- gsub("XX", "AS", colnames(met) )
counts <- counts[,!(names(counts) %in% c("VS16_6AS_D2","Exp2_F1_7AS_D3","WBC3_5AS_D1","Exp1_F1_12AS_D3","DAC6_6AS_D1","WBC3_8AS_D2","VS16_6AS_D3"))]
counts <- counts[,colnames(counts) %in% colnames(met)]
met <- met[,colnames(met) %in% colnames(counts)]
write.csv(met,'Metabolites_matrix_XX_exp1_overlapping_with_AS_expression.csv',row.names=T,quote=F)

# VS16
VS16 = counts[,grep(colnames(counts),pattern='VS16',fixed = TRUE)]
group <- factor(substr(colnames(VS16),nchar(colnames(VS16))-1,nchar(colnames(VS16))))
d <- DGEList(counts=VS16,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("D1", "D2"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('VS16_D2_vs_D1_',colnames(res1),sep='')
out <- res1

# WBC3
WBC3 = counts[,grep(colnames(counts),pattern='WBC3',fixed = TRUE)]
group <- factor(substr(colnames(WBC3),nchar(colnames(WBC3))-1,nchar(colnames(WBC3))))
d <- DGEList(counts=WBC3,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("D1", "D2"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('WBC3_D2_vs_D1_',colnames(res1),sep='')
out <- cbind(out[order(rownames(out)),],res1[order(rownames(res1)),])

# Exp1_F1
Exp1_F1 = counts[,grep(colnames(counts),pattern='Exp1_F1',fixed = TRUE)]
group <- factor(substr(colnames(Exp1_F1),nchar(colnames(Exp1_F1))-1,nchar(colnames(Exp1_F1))))
d <- DGEList(counts=Exp1_F1,group=group)
d <- calcNormFactors(d)
d$samples
levels(d$samples$group)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
FG_de.com <- exactTest(d, pair=c("D1", "D2"))
results_FG <- topTags(FG_de.com,n = length(d$AveLogCPM))
res1 <- results_FG$table[,c(1,4)]
colnames(res1) <- paste('Exp1_F1_D2_vs_D1_',colnames(res1),sep='')
out <- cbind(out[order(rownames(out)),],res1[order(rownames(res1)),])

write.csv(out,'FC_Exp1_linking_to_XX_met.csv',row.names=T,quote=F)
FC <- out[(abs(out[,1]) > 1 & out[,2] < 0.05) | (abs(out[,3]) > 1 & out[,4] < 0.05) | (abs(out[,5]) > 1 & out[,6] < 0.05), ]
write.csv(FC,'FC_Exp1_linking_to_XX_met_responsive_in_at_least_one_sample.csv',row.names=T,quote=F)



#######################################################################################
###############															###############
###############			Get median TPM and metabolic   					###############
###############				                       		   				###############
#######################################################################################
setwd('D:\\Projects\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcirptome_Xingxing')
TPM = read.table("All_but_AP13_TPM_unique_AS.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
TPM <- TPM[,grep(colnames(TPM),pattern='AS',fixed = TRUE)]

# get median values across replicates
library(matrixStats)
met_2 <- read.csv('Metabolites_matrix_XX_exp1_overlapping_with_AS_expression.csv',row.names=1,head=T,stringsAsFactors=F)
#TPM_2 <- read.table('Exp_matrix_XX_exp1.txt',row.names=1,head=T,stringsAsFactors=F)
TPM_2 <- TPM[,colnames(TPM) %in% colnames(met_2)]
# TPM_median <- c()
# for(cul in c('VS16','WBC3','F1')){
	# subdat <- TPM_2[,grep(colnames(TPM_2),pattern=cul,fixed = TRUE)]
	# for(days in c('D1','D2')){
		# subtem <- subdat[,grep(colnames(subdat),pattern=days,fixed = TRUE)]
		# TPM_median <- cbind(TPM_median,rowMedians(as.matrix(subtem)))
		# colnames(TPM_median)[dim(TPM_median)[2]] = paste(cul,days,sep='_')
		# }
	# }

met_median <- c()
for(cul in c('VS16','WBC3','F1')){
	subdat <- met_2[,grep(colnames(met_2),pattern=cul,fixed = TRUE)]
	for(days in c('D1','D2')){
		subtem <- subdat[,grep(colnames(subdat),pattern=days,fixed = TRUE)]
		met_median <- cbind(met_median,rowMedians(as.matrix(subtem)))
		colnames(met_median)[dim(met_median)[2]] = paste(cul,days,sep='_')
		}
	}
rownames(met_median) <- rownames(met_2)
write.table(met_median,'Metabolites_matrix_XX_exp1_overlapping_with_AS_expression_median.txt',row.names=T,quote=F,sep='\t')
#write.table(TPM_median,'TPM_AS_Exp1_overlapping_with_XX_median.txt',row.names=T,quote=F,sep='\t')

########################################################################################################
###########################################################
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
# filter genes with small variation across samples. 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
#m.mad <- apply(dataExpr,1,mad)
#dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
# soft threshold
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=70, by=2))
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('TPM_xingxing_Exp1_expressed_differential_median_soft.pdf')
par(mfrow = c(2,1))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
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
# 筛选标准。R-square=0.9
dev.off()

# build network
enableWGCNAThreads()
softPower = 32
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
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf('Correlation_among_AS_exp1_overlapping_with_XX_met_WGCNA_modules_median_MEDissThres_0.1.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
#########################
# associate the expression modules with traits
trait <- "Metabolites_matrix_XX_exp1_overlapping_with_AS_expression_median.txt"
# 读入表型数据，不是必须的
traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
colnames(traitData) <- rownames(dataExpr)						
traitData <- t(traitData)

### 模块与表型数据关联
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

# signif表示保留几位小数
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
# filter genes with small variation across samples. 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
#m.mad <- apply(dataExpr,1,mad)
#dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
# soft threshold
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=70, by=2))
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('Metabolites_XX_Exp1_overlapping_with_AS_expression_soft.pdf')
par(mfrow = c(2,1))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
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
# 筛选标准。R-square=0.85
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
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf('Correlation_among_AS_exp1_overlapping_with_XX_met_WGCNA_modules_median_MEDissThres_0.1.pdf')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

# zscore normalization
# https://www.rdocumentation.org/packages/rcompanion/versions/2.4.6/topics/blom
# remove genes with all TPM < 2





















# res <- c()
# for(i in 1:nrow(TPM_2)){
	# for(j in 1:nrow(met_2)){
		# res <- rbind(res,c(rownames(TPM_2)[i],rownames(met_2)[j],cor.test(TPM_2[i,],met_2[j,])[4]))
		# }
	# if(i%%100 == 0) print(i)
	# }
res <- cor(t(TPM_2),t(met_2))
res2 <- na.omit(res)
res2 <- round(res2,3)
write.csv(res2,'Correlation_between_metabolome_transscriptome_exp1_XX.csv',row.names=TRUE,quote=F)

TPM_3 <- TPM_2[apply(TPM_2 >= 2, 1, any),]
res3 <- cor(t(TPM_3),t(met_2))
res3 <- na.omit(res3)
res4 <- round(res3,3)
write.csv(res4,'Correlation_between_metabolome_transscriptome_exp1_XX_TPM_2.csv',row.names=TRUE,quote=F)
tem <- res4[order(res4[,1],decreasing=T),]

# corre <- res4[apply(abs(res4) >= 0.5, 1, any),]

# res5 <- cor(t(TPM_3),t(met_2),method='spearman')
# res5 <- na.omit(res5)
# res6 <- round(res5,3)
# write.csv(res6,'Correlation_between_metabolome_transscriptome_exp1_XX_TPM_2_Spearman.csv',row.names=F,quote=F)
# tem2 <- res6[order(res6[,1],decreasing=T),]

# met_2 <- cbind(met_2,rowSums(met_2))
# met_2 <- met_2[order(met_2[,ncol(met_2)],decreasing=T),]
# tem2 <- res6[order(res6[,2208],decreasing=T),]

# tem <- res4[order(res4[,2208],decreasing=T),]


met_2_D1 <- met_2[,grep(colnames(met_2),pattern='D1',fixed = TRUE)]
met_2_D2 <- met_2[,grep(colnames(met_2),pattern='D2',fixed = TRUE)]
TPM_3 <- TPM_2[apply(TPM_2 >= 2, 1, any),]
TPM_3_D1 <- TPM_3[,grep(colnames(TPM_3),pattern='D1',fixed = TRUE)]
TPM_3_D2 <- TPM_3[,grep(colnames(TPM_3),pattern='D2',fixed = TRUE)]
# D1
group <- factor(substr(colnames(met_2_D1),1,4))
out <- c()
for(i in 1:nrow(met_2_D1)){
	out <- rbind(out,c(rownames(met_2_D1)[i],kruskal.test(met_2_D1[i,], group)[3]))
	}
out <- out[order(as.numeric(out[,2])),]
out <- cbind(out,p.adjust(as.numeric(out[,2]), method = "BH"))
subout <- as.data.frame(out[as.numeric(out[,3]) < 0.05,]) # 1118
write.csv(subout,'Significant_diff_metabolites.csv',row.names=TRUE,quote=F)

correl_D1 <- cor(t(TPM_3_D1),t(met_2_D1[as.character(subout$V1),]))
write.csv(correl_D1,'Correlation_between_sig_dif_metabolome_transscriptome_exp1_XX_TPM_2_D1.csv',row.names=TRUE,quote=F)

########################################
# D2 vs D1, significantly different abundance for metabolites
# D2 vs D1, vs16
met_2_vs16 <- met_2[,grep(colnames(met_2),pattern='VS16',fixed = TRUE)]
group <- factor(substr(colnames(met_2_vs16),nchar(colnames(met_2_vs16))-1,nchar(colnames(met_2_vs16))))
out <- c()
for(i in 1:nrow(met_2_vs16)){
	#out <- rbind(out,c(rownames(met_2_vs16)[i],kruskal.test(met_2_vs16[i,], group)[3]))
	out <- rbind(out,c(rownames(met_2_vs16)[i],wilcox.test(met_2_vs16[i,1:10],met_2_vs16[i,11:19])[3]))
	}
out <- out[order(as.numeric(out[,2])),]
out <- cbind(out,p.adjust(as.numeric(out[,2]), method = "BH"))
subout <- out[as.numeric(out[,3]) < 0.05,] # 10
write.csv(subout,'Sig_diff_metabolites_VS16_D2_vs_D1.csv',quote=F,row.names=F)


# D2 vs D1, wbc3
met_2_WBC3 <- met_2[,grep(colnames(met_2),pattern='WBC3',fixed = TRUE)]
group <- factor(substr(colnames(met_2_WBC3),nchar(colnames(met_2_WBC3))-1,nchar(colnames(met_2_WBC3))))
out <- c()
for(i in 1:nrow(met_2_WBC3)){
	#out <- rbind(out,c(rownames(met_2_WBC3)[i],kruskal.test(met_2_WBC3[i,], group)[3]))
	out <- rbind(out,c(rownames(met_2_WBC3)[i],wilcox.test(met_2_WBC3[i,1:10],met_2_WBC3[i,11:20])[3]))
	}
out <- out[order(as.numeric(out[,2])),]
out <- cbind(out,p.adjust(as.numeric(out[,2]), method = "BH"))
subout <- out[as.numeric(out[,3]) < 0.05,] # 12
write.csv(subout,'Sig_diff_metabolites_WBC3_D2_vs_D1.csv',quote=F,row.names=F)

# D2 vs D1, F1
met_2_F1 <- met_2[,grep(colnames(met_2),pattern='F1',fixed = TRUE)]
group <- factor(substr(colnames(met_2_F1),nchar(colnames(met_2_F1))-1,nchar(colnames(met_2_F1))))
out <- c()
for(i in 1:nrow(met_2_F1)){
	#out <- rbind(out,c(rownames(met_2_F1)[i],kruskal.test(met_2_F1[i,], group)[3]))
	out <- rbind(out,c(rownames(met_2_F1)[i],wilcox.test(met_2_F1[i,1:14],met_2_F1[i,15:28])[3]))
	}
out <- out[order(as.numeric(out[,2])),]
out <- cbind(out,p.adjust(as.numeric(out[,2]), method = "BH"))
subout <- out[as.numeric(out[,3]) < 0.05,] # 62
write.csv(subout,'Sig_diff_metabolites_Exp1_F1_D2_vs_D1.csv',quote=F,row.names=F)

######
setwd('D:\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcriptome')
dat <- read.csv('Sig_diff_metabolites_XX_Exp1_D2_vs_D1.csv')
F1 <- dat[dat[,4]=='F1',]
WBC3 <- dat[dat[,4]=='WBC3',]
VS16 <- dat[dat[,4]=='VS16',]
library(gplots)
pdf("Sig_diff_metabolites_XX_Exp1_D2_vs_D1.pdf")
venn(list(VS16[,1],WBC3[,1],F1[,1]))
dev.off()
shared <- intersect(F1[,1],intersect(VS16[,1],WBC3[,1]))

TPM = read.table("All_but_AP13_TPM_unique.txt", sep="\t", row.names=1,head=T,stringsAsFactors=F)
TPM <- TPM[,grep(colnames(TPM),pattern='XX',fixed = TRUE)]
TPM_t <- t(TPM)
ID <- read.csv("Exp1_metabolites_XX_sample_ID.csv")
met <- read.csv("Exp1_metabolites_XX.csv",row.names=1)
met_t <- t(met)
met_t_rowname <- rownames(met_t)
for(i in 1:nrow(met_t)){
	subdat <- ID[ID[,1]==met_t_rowname[i],]
	rownames(met_t)[i] = as.character(subdat[1,3])
	}

TPM_t <- TPM_t[intersect(rownames(met_t),rownames(TPM_t)),]
met_t <- met_t[intersect(rownames(met_t),rownames(TPM_t)),]
TPM_2 <- t(TPM_t)
met_2 <- t(met_t)

met_F1 <- met_2[F1[,1],]
met_F1 <- log2(met_F1*100+1)
library(gplots)
pdf('metabolic_sig_different_F1_Exp1.pdf')
heatmap.2(as.matrix(met_F1),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

met_WBC3 <- met_2[WBC3[,1],]
met_WBC3 <- log2(met_WBC3*100+1)
library(gplots)
pdf('metabolic_sig_different_WBC3_Exp1.pdf')
heatmap.2(as.matrix(met_WBC3),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

met_VS16 <- met_2[VS16[,1],]
met_VS16 <- log2(met_VS16*100+1)
library(gplots)
pdf('metabolic_sig_different_VS16_Exp1.pdf')
heatmap.2(as.matrix(met_VS16),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

met_vs16 <- setdiff(VS16[,1], union(F1[,1],WBC3[,1]))
met_VS16_specific <- met_2[met_vs16,]
met_VS16_specific <- log2(met_VS16_specific*100+1)
library(gplots)
pdf('metabolic_sig_different_VS16_specific_Exp1.pdf')
heatmap.2(as.matrix(met_VS16_specific),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

met_all <- met_2[unique(dat[,1]),]
met_all <- log2(met_all*100+1)
library(gplots)
pdf('metabolic_sig_different_all_Exp1.pdf')
heatmap.2(as.matrix(met_all),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()



###### correlation between metabolome and transcriptome for only the 3 shared metabolites among VS16, WBC3, and F1
rownames(TPM_2) <- gsub("Pavir.","",rownames(TPM_2))
rownames(TPM_2) <- gsub(".v5.1","",rownames(TPM_2))
TPM_10 <- TPM_2[apply(TPM_2 >= 2, 1, any),]
met_2 <- t(met_t)
met_2_common <- met_2[c("13.76_256.2635m/z", "13.99_281.2718n", "14.01_264.2453n"),]
###### pcc
pcc <- cor(t(TPM_10),t(met_2_common))
pcc <- na.omit(pcc)
pcc <- round(pcc,3)
write.csv(pcc,'Correlation_between_metabolome_transscriptome_exp1_XX_TPM_2_3_common_sig_dif_metabolites.csv',row.names=TRUE,quote=F)
pcc <- pcc[order(pcc[,3],decreasing=T),]
pcc_0.8 <- pcc[pcc[,3] > 0.8,]
TPM_all <- TPM_2[rownames(pcc_0.8),]
TPM_all <- log2(TPM_all + 1)
pdf('TPM_for_genes_correlated_with_metabolic_sig_different_all_Exp1_pcc0.8.pdf')
heatmap.2(as.matrix(TPM_all),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = FALSE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

# draw the heatmap of the TPM
pcc_0.7 <- pcc[pcc[,3] > 0.7,]
TPM_all <- TPM_2[rownames(pcc_0.7),]
TPM_all <- log2(TPM_all + 1)
pdf('TPM_for_genes_correlated_with_metabolic_sig_different_all_Exp1_pcc0.7.pdf')
heatmap.2(as.matrix(TPM_all),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = FALSE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
write.table(rownames(TPM_all),'Genes_correlated_with_metabolic_sig_different_all_Exp1_pcc0.7.txt',row.names=F,col.names=F,quote=F)
GO_pipeline('Genes_correlated_with_metabolic_sig_different_all_Exp1_pcc0.7.txt','Genes_correlated_with_metabolic_sig_different_all_Exp1_pcc0.7')


###### spearman
pcc <- cor(t(TPM_10),t(met_2_common),method='spearman')
pcc <- na.omit(pcc)
pcc <- round(pcc,3)
write.csv(pcc,'Correlation_between_metabolome_transscriptome_exp1_XX_TPM_2_3_common_sig_dif_metabolites_spearman.csv',row.names=TRUE,quote=F)
pcc <- pcc[order(pcc[,3],decreasing=T),]
pcc_0.7 <- pcc[pcc[,3] > 0.7,]
TPM_all <- TPM_2[rownames(pcc_0.7),]
TPM_all <- log2(TPM_all + 1)
pdf('TPM_for_genes_correlated_with_metabolic_sig_different_all_Exp1_spearman0.7.pdf')
heatmap.2(as.matrix(TPM_all),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = FALSE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
write.table(rownames(TPM_all),'Genes_correlated_with_metabolic_sig_different_all_Exp1_spearman0.7.txt',row.names=F,col.names=F,quote=F)
GO_pipeline('Genes_correlated_with_metabolic_sig_different_all_Exp1_spearman0.7.txt','Genes_correlated_with_metabolic_sig_different_all_Exp1_spearman0.7')

#################################################################################################
################                                                                 ################
################                       RNA seq analysis                          ################
################                                                                 ################
#################################################################################################
############################################################################################
############### GO enrichment for differentially expressed genes
setwd('D:\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcriptome')
FC <- read.table('D:\\Switchgrass_allelic_specific\\03_RNA_seq\\All_switchgrass_cold_FC_unique_AS.txt')
TPM <- read.table('D:\\Switchgrass_allelic_specific\\03_RNA_seq\\All_switchgrass_cold_median_TPM_unique_AS.txt')
# expressed <- TPM[apply(TPM > 2,1,any),]
# FC <- FC[rownames(expressed),]
rownames(FC) <- gsub("Pavir.","",rownames(FC))
rownames(FC) <- gsub(".v5.1","",rownames(FC))
VS16_up <- FC[FC$VS16_AS_D2_vs_D1_logFC > 1 & FC$VS16_AS_D2_vs_D1_FDR < 0.05,]
write.table(rownames(VS16_up),'VS16_AS_D2_vs_D1_up_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

WBC3_up <- FC[FC$WBC3_AS_D2_vs_D1_logFC > 1 & FC$WBC3_AS_D2_vs_D1_FDR < 0.05,]
write.table(rownames(WBC3_up),'WBC3_AS_D2_vs_D1_up_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

AP13_up <- FC[FC$AP13_D2_vs_D1_logFC > 1 & FC$AP13_D2_vs_D1_FDR < 0.05,]
write.table(rownames(AP13_up),'AP13_AS_D2_vs_D1_up_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

DAC6_up <- FC[FC$DAC6_D2_vs_D1_logFC > 1 & FC$DAC6_D2_vs_D1_FDR < 0.05,]
write.table(rownames(DAC6_up),'DAC6_AS_D2_vs_D1_up_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

VS16_down <- FC[FC$VS16_AS_D2_vs_D1_logFC < -1 & FC$VS16_AS_D2_vs_D1_FDR < 0.05,]
write.table(rownames(VS16_down),'VS16_AS_D2_vs_D1_down_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

WBC3_down <- FC[FC$WBC3_AS_D2_vs_D1_logFC < -1 & FC$WBC3_AS_D2_vs_D1_FDR < 0.05,]
write.table(rownames(WBC3_down),'WBC3_AS_D2_vs_D1_down_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

AP13_down <- FC[FC$AP13_D2_vs_D1_logFC < -1 & FC$AP13_D2_vs_D1_FDR < 0.05,]
write.table(rownames(AP13_down),'AP13_AS_D2_vs_D1_down_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')

DAC6_down <- FC[FC$DAC6_D2_vs_D1_logFC < -1 & FC$DAC6_D2_vs_D1_FDR < 0.05,]
write.table(rownames(DAC6_down),'DAC6_AS_D2_vs_D1_down_regulated_genes.txt',row.names=F,col.names=F,quote=F,sep='\t')



Enrichment <- function(k,n,C,G){
	return((k / C) / (n/ G))
	}
library('GOstats')
library('GSEABase')
library(pvclust)
library(gplots)
GO <- read.table('D:\\Transformer_switchgrass\\05_GO_enrichment\\Switchgrass_GO_annotation.txt',head=F,sep='\t')
# null <- rownames(read.table('D:\\Transformer_switchgrass\\02_RNA-seq_results\\Swithchgrass_genes_FC_nochange_FC0.5_FDR0.1_unique.txt',head=T,sep='\t'))
# null <- gsub("Pavir.","",null)
# null <- gsub(".v5.1","",null)
# null_go <- unique(merge(GO,t(t(null)),by.x='V1',by.y='V1'))
# null_num <- length(null)
null_go <- unique(GO)
null_num <- 80278

GO_pipeline <- function(Genes_file,short_name_to_save){
	pos <- read.table(Genes_file,head=F,sep='\t')
	pos[,1] <- gsub("Pavir.","",pos[,1])
	pos[,1] <- gsub(".v5.1","",pos[,1])
	pos_go <- unique(merge(GO,pos,by.x='V1',by.y='V1'))
	pos_num <- nrow(pos)
	res <- c()
	go <- unique(GO[,2])
	for(i in 1:length(go)){
		Pos <- pos_go[pos_go[,2]==go[i],]
		Null <- null_go[null_go[,2]==go[i],]
		a = nrow(Pos) # genes in cluster and with GO term
		b = pos_num - a # genes in cluster but with no GO term
		cc = nrow(Null) # genes with GO but not in cluster
		d = null_num - cc# genes with no GO and also not in cluster
		out <- c(go[i],a,b,cc,d)
		if(a+cc != 0) res <- rbind(res,out)
		}
	colnames(res) <- c('GO','Genes_responsive_with_GO','Genes_responsive_with_no_GO','Genes_not_responsive_with_GO','Genes_not_responsive_with_no_GO')
	write.table(res,paste('GO_',short_name_to_save,'.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	enrichment <- c()
	for(i in 1:nrow(res)){
		numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
		p <- fisher.test(numbers, alternative = "two.sided")[[1]][1]
		a = as.numeric(res[i,2])
		b = as.numeric(res[i,3])
		cc = as.numeric(res[i,4])
		d = as.numeric(res[i,5])
		if(Enrichment(a, a+b, a+cc, a+b+cc+d) >= 1) direction = '+' else direction = '-'
		enrichment <- rbind(enrichment, c(res[i,],direction,p))
		}
	write.table(enrichment,paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- read.table(paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
	dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
	dat <- dat[order(dat[,8]),]  
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- cbind(dat,'BP'='')
	dat <- cbind(dat,'CC'='')
	dat <- cbind(dat,'MF'='')
	for(i in 1:nrow(dat)){
		tryCatch( {
			if(!is.null(getGOTerm(dat[i,1])$BP[1])) dat[i,9] <- getGOTerm(dat[i,1])$BP[1]
			if(!is.null(getGOTerm(dat[i,1])$CC[1])) dat[i,10] <- getGOTerm(dat[i,1])$CC[1]
			if(!is.null(getGOTerm(dat[i,1])$MF[1])) dat[i,11] <- getGOTerm(dat[i,1])$MF[1]
			},
			error = function(e) {print(paste("no GO for ",dat[i,1]));NaN},
			finally = {})
		}
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	
	subdat <- dat[dat[,8] < 0.05 & dat$BP != '',]
	subdat <- cbind(subdat,'logP'=0)
	for(i in 1:nrow(subdat)){
		if(subdat[i,6] == '-')subdat[i,ncol(subdat)] = log10(subdat[i,8])
		if(subdat[i,6] == '+')subdat[i,ncol(subdat)] = -log10(subdat[i,8])
		}
	subdat[subdat$logP < -10,12] <- -10
	subdat[subdat$logP > 10,12] <- 10
	df <- as.matrix(cbind(subdat$logP,c(10,-10)))
	rownames(df) <- subdat$BP
	df <- df[order(df[,1],decreasing=T),]
	pdf(paste("GO_", short_name_to_save,'.fisher.qvalue_GO_term_BP.pdf',sep=''))
	heatmap.2(df,
		trace="none", # remove trace lines
		col = colorRampPalette(c("blue","white","red"))(21), # color palette
		dendrogram="none", # remove dendrogram
		cexRow=0.8,# row label font size
		labCol=F, # remove column labels
		Rowv=F, # turn off row clustering
		Colv=F, # turn off column clustering
		lmat=rbind(4:3,2:1),
		lwid=c(0.8,4), # dimensions of display array cell widths
		lhei=c(0.8,4)) # dimensions of display cell heights
	dev.off()
	write.table(subdat[,c(9,8,ncol(subdat))],paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term_BP.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')

	}


GO_pipeline('VS16_AS_D2_vs_D1_up_regulated_genes.txt','VS16_AS_D2_vs_D1_up')
GO_pipeline('WBC3_AS_D2_vs_D1_up_regulated_genes.txt','WBC3_AS_D2_vs_D1_up')
GO_pipeline('AP13_AS_D2_vs_D1_up_regulated_genes.txt','AP13_AS_D2_vs_D1_up')
GO_pipeline('DAC6_AS_D2_vs_D1_up_regulated_genes.txt','DAC6_AS_D2_vs_D1_up')

GO_pipeline('VS16_AS_D2_vs_D1_down_regulated_genes.txt','VS16_AS_D2_vs_D1_down')
GO_pipeline('WBC3_AS_D2_vs_D1_down_regulated_genes.txt','WBC3_AS_D2_vs_D1_down')
GO_pipeline('AP13_AS_D2_vs_D1_down_regulated_genes.txt','AP13_AS_D2_vs_D1_down')
GO_pipeline('DAC6_AS_D2_vs_D1_down_regulated_genes.txt','DAC6_AS_D2_vs_D1_down')

VS16 <- read.table('GO_VS16_AS_D2_vs_D1_up.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
WBC3 <- read.table('GO_WBC3_AS_D2_vs_D1_up.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
AP13 <- read.table('GO_AP13_AS_D2_vs_D1_up.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
DAC6 <- read.table('GO_DAC6_AS_D2_vs_D1_up.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')

library(gplots)
pdf("GO_BP_upregulated_venn.pdf")
venn(list(VS16[,1],WBC3[,1],AP13[,1],DAC6[,1]))
dev.off()

intersect(VS16[,1],WBC3[,1])
intersect(AP13[,1],DAC6[,1])
intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1]))

setdiff(intersect(VS16[,1],DAC6[,1]),union(WBC3[,1],AP13[,1]))
setdiff(intersect(AP13[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(intersect(WBC3[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(WBC3[,1],union(union(AP13[,1],DAC6[,1]),VS16[,1]))
setdiff(AP13[,1],union(union(WBC3[,1],DAC6[,1]),VS16[,1]))
setdiff(intersect(AP13[,1],WBC3[,1]),union(DAC6[,1],VS16[,1]))
setdiff(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(DAC6[,1],union(union(AP13[,1],WBC3[,1]),VS16[,1]))

intersect(intersect(VS16[,1],WBC3[,1]),AP13[,1])
intersect(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(intersect(DAC6[,1],AP13[,1]),WBC3[,1])


VS16 <- read.table('GO_VS16_AS_D2_vs_D1_down.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
WBC3 <- read.table('GO_WBC3_AS_D2_vs_D1_down.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
AP13 <- read.table('GO_AP13_AS_D2_vs_D1_down.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')
DAC6 <- read.table('GO_DAC6_AS_D2_vs_D1_down.fisher.qvalue_GO_term_BP.txt',head=T,sep='\t')

library(gplots)
pdf("GO_BP_downregulated_venn.pdf")
venn(list(VS16[,1],WBC3[,1],AP13[,1],DAC6[,1]))
dev.off()

intersect(VS16[,1],WBC3[,1])
intersect(AP13[,1],DAC6[,1])
intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1]))

setdiff(intersect(VS16[,1],DAC6[,1]),union(WBC3[,1],AP13[,1]))
setdiff(intersect(AP13[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(intersect(WBC3[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(WBC3[,1],union(union(AP13[,1],DAC6[,1]),VS16[,1]))
setdiff(AP13[,1],union(union(WBC3[,1],DAC6[,1]),VS16[,1]))
setdiff(intersect(AP13[,1],WBC3[,1]),union(DAC6[,1],VS16[,1]))
setdiff(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(DAC6[,1],union(union(AP13[,1],WBC3[,1]),VS16[,1]))

intersect(intersect(VS16[,1],WBC3[,1]),AP13[,1])
intersect(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(intersect(DAC6[,1],AP13[,1]),WBC3[,1])


PATHWAY <- read.table('D:\\Transformer_switchgrass\\08_pathway_enrichment\\Switchgrass_pathway_annotation.txt',head=F,sep='\t')
null_pathway <- unique(PATHWAY)
null_num <- 80278
Pathway_pipeline <- function(Genes_file,short_name_to_save){
	pos <- read.table(Genes_file,head=F,sep='\t')
	pos[,1] <- gsub("Pavir.","",pos[,1])
	pos[,1] <- gsub(".v5.1","",pos[,1])
	pos_pathway <- unique(merge(PATHWAY,pos,by.x='V1',by.y='V1'))
	pos_num <- nrow(pos)
	res <- c()
	pathway <- unique(PATHWAY[,2])
	for(i in 1:length(pathway)){
		Pos <- pos_pathway[pos_pathway[,2]==pathway[i],]
		Null <- null_pathway[null_pathway[,2]==pathway[i],]
		a = nrow(Pos) # genes in cluster and with pathway term
		b = pos_num - a # genes in cluster but with no pathway term
		cc = nrow(Null) # genes with pathway but not in cluster
		d = null_num - cc# genes with no pathway and also not in cluster
		out <- c(pathway[i],a,b,cc,d)
		if(a+cc != 0) res <- rbind(res,out)
		}
	colnames(res) <- c('Pathway','Genes_responsive_in_pathway','Genes_responsive_not_in_pathway','Genes_not_responsive_in_pathway','Genes_not_responsive_not_in_pathway')
	write.table(res,paste('Pathway_',short_name_to_save,'.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	enrichment <- c()
	for(i in 1:nrow(res)){
		numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
		p <- fisher.test(numbers, alternative = "two.sided")[[1]][1]
		a = as.numeric(res[i,2])
		b = as.numeric(res[i,3])
		cc = as.numeric(res[i,4])
		d = as.numeric(res[i,5])
		if(Enrichment(a, a+b, a+cc, a+b+cc+d) >= 1) direction = '+' else direction = '-'
		enrichment <- rbind(enrichment, c(res[i,],direction,p))
		}
	write.table(enrichment,paste('Pathway_',short_name_to_save,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- read.table(paste('Pathway_',short_name_to_save,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
	dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
	dat <- dat[order(dat[,8]),]  
	write.table(dat,paste('Pathway_',short_name_to_save,'.fisher.qvalue.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')
	
	subdat <- dat[dat[,8] < 0.05,]
	subdat <- cbind(subdat,'logP'=0)
	for(i in 1:nrow(subdat)){
		if(subdat[i,6] == '-')subdat[i,ncol(subdat)] = log10(subdat[i,8])
		if(subdat[i,6] == '+')subdat[i,ncol(subdat)] = -log10(subdat[i,8])
		}
	subdat[subdat$logP < -10,9] <- -10
	subdat[subdat$logP > 10,9] <- 10
	df <- as.matrix(cbind(subdat$logP,c(10,-10)))
	rownames(df) <- subdat$V1
	df <- df[order(df[,1],decreasing=T),]
	pdf(paste("Pathway_", short_name_to_save,'.fisher.qvalue.pdf',sep=''))
	heatmap.2(df,
		trace="none", # remove trace lines
		col = colorRampPalette(c("blue","white","red"))(21), # color palette
		dendrogram="none", # remove dendrogram
		cexRow=0.8,# row label font size
		labCol=F, # remove column labels
		Rowv=F, # turn off row clustering
		Colv=F, # turn off column clustering
		lmat=rbind(4:3,2:1),
		lwid=c(0.8,4), # dimensions of display array cell widths
		lhei=c(0.8,4)) # dimensions of display cell heights
	dev.off()
	write.table(subdat[,c(1,8,9)],paste('Pathway_',short_name_to_save,'.fisher.qvalue.txt',sep=''),row.names=T,col.names=T,quote=F,sep='\t')
	}

Pathway_pipeline('VS16_AS_D2_vs_D1_up_regulated_genes.txt','VS16_AS_D2_vs_D1_up')
Pathway_pipeline('WBC3_AS_D2_vs_D1_up_regulated_genes.txt','WBC3_AS_D2_vs_D1_up')
Pathway_pipeline('AP13_AS_D2_vs_D1_up_regulated_genes.txt','AP13_AS_D2_vs_D1_up')
Pathway_pipeline('DAC6_AS_D2_vs_D1_up_regulated_genes.txt','DAC6_AS_D2_vs_D1_up')

VS16 <- read.table('Pathway_VS16_AS_D2_vs_D1_up.fisher.qvalue.txt',head=T,sep='\t')
WBC3 <- read.table('Pathway_WBC3_AS_D2_vs_D1_up.fisher.qvalue.txt',head=T,sep='\t')
AP13 <- read.table('Pathway_AP13_AS_D2_vs_D1_up.fisher.qvalue.txt',head=T,sep='\t')
DAC6 <- read.table('Pathway_DAC6_AS_D2_vs_D1_up.fisher.qvalue.txt',head=T,sep='\t')

library(gplots)
pdf("Pathway_upregulated_venn.pdf")
venn(list(VS16[,1],WBC3[,1],AP13[,1],DAC6[,1]))
dev.off()

intersect(VS16[,1],WBC3[,1])
intersect(AP13[,1],DAC6[,1])
intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1]))

setdiff(intersect(VS16[,1],DAC6[,1]),union(WBC3[,1],AP13[,1]))
setdiff(intersect(AP13[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(intersect(WBC3[,1],intersect(VS16[,1],DAC6[,1])),intersect(intersect(VS16[,1],WBC3[,1]),intersect(AP13[,1],DAC6[,1])))
setdiff(WBC3[,1],union(union(AP13[,1],DAC6[,1]),VS16[,1]))
setdiff(AP13[,1],union(union(WBC3[,1],DAC6[,1]),VS16[,1]))
setdiff(intersect(AP13[,1],WBC3[,1]),union(DAC6[,1],VS16[,1]))
setdiff(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(DAC6[,1],union(union(AP13[,1],WBC3[,1]),VS16[,1]))

intersect(intersect(VS16[,1],WBC3[,1]),AP13[,1])
intersect(intersect(DAC6[,1],WBC3[,1]),AP13[,1])
setdiff(intersect(DAC6[,1],AP13[,1]),WBC3[,1])


###############################################################
# check the expression of genes with enriched GOs
library(gplots)
TPM <- read.table('D:\\Switchgrass_allelic_specific\\03_RNA_seq\\All_switchgrass_cold_median_TPM_unique_AS.txt')
rownames(TPM) <- gsub("Pavir.","",rownames(TPM))
rownames(TPM) <- gsub(".v5.1","",rownames(TPM))
GO <- read.table('D:\\Transformer_switchgrass\\05_GO_enrichment\\Switchgrass_GO_annotation.txt',head=F,sep='\t')
GO_target <- GO[GO[,2]=='GO:0005992',]
TPM_target <- TPM[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('trehalose_biosynthetic_process_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

GO_target <- GO[GO[,2]=='GO:0001522',]
TPM_target <- TPM[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('pseudouridine_synthesis_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()

GO_target <- GO[GO[,2]=='GO:0007165',]
TPM_target <- TPM[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('signal_transduction_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
VS16_target <- merge(VS16,GO_target,by.x='V1',by.y='V1')
WBC3_target <- merge(WBC3,GO_target,by.x='V1',by.y='V1')
AP13_target <- merge(AP13,GO_target,by.x='V1',by.y='V1')
DAC6_target <- merge(DAC6,GO_target,by.x='V1',by.y='V1')

###############################################################
# venn digram of up-regulated and down-regulated genes
VS16 <- read.table('VS16_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
WBC3 <- read.table('WBC3_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
AP13 <- read.table('AP13_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
DAC6 <- read.table('DAC6_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
pdf("Up-regulated_genes_venn.pdf")
venn(list(VS16[,1],WBC3[,1],AP13[,1],DAC6[,1]))
dev.off()

upland <- setdiff(intersect(VS16[,1],DAC6[,1]),union(WBC3[,1],AP13[,1]))
write.table(upland,'Genes_upregulated_only_in_upland.txt',row.names=F,col.names=F,quote=F,sep='\t')
GO_pipeline('Genes_upregulated_only_in_upland.txt','Genes_upregulated_only_in_upland')

VS16 <- read.table('VS16_AS_D2_vs_D1_down_regulated_genes.txt',head=F,sep='\t')
WBC3 <- read.table('WBC3_AS_D2_vs_D1_down_regulated_genes.txt',head=F,sep='\t')
AP13 <- read.table('AP13_AS_D2_vs_D1_down_regulated_genes.txt',head=F,sep='\t')
DAC6 <- read.table('DAC6_AS_D2_vs_D1_down_regulated_genes.txt',head=F,sep='\t')
pdf("Down-regulated_genes_venn.pdf")
venn(list(VS16[,1],WBC3[,1],AP13[,1],DAC6[,1]))
dev.off()
###


####################
# all up-regulated genes in D2
VS16 <- read.table('VS16_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
WBC3 <- read.table('WBC3_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
AP13 <- read.table('AP13_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
DAC6 <- read.table('DAC6_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
all_genes <- union(union(VS16[,1],WBC3[,1]),union(AP13[,1],DAC6[,1]))
TPM_expressed <- TPM[all_genes,]
GO <- GO[GO[,1] %in% all_genes,]
GO_target <- GO[GO[,2]=='GO:0005992',]
TPM_target <- TPM_expressed[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('trehalose_biosynthetic_process_up-regulated_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
VS16_target <- merge(VS16,GO_target,by.x='V1',by.y='V1')
WBC3_target <- merge(WBC3,GO_target,by.x='V1',by.y='V1')
AP13_target <- merge(AP13,GO_target,by.x='V1',by.y='V1')
DAC6_target <- merge(DAC6,GO_target,by.x='V1',by.y='V1')


GO_target <- GO[GO[,2]=='GO:0001522',]
TPM_target <- TPM_expressed[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('pseudouridine_synthesis_up-regulated_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()


GO_target <- GO[GO[,2]=='GO:0007165',]
TPM_target <- TPM_expressed[GO_target[,1],]
TPM_target <- log2(TPM_target+1)
pdf('signal_transduction_up-regulated_genes_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()




VS16 <- read.table('VS16_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
WBC3 <- read.table('WBC3_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
AP13 <- read.table('AP13_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
DAC6 <- read.table('DAC6_AS_D2_vs_D1_up_regulated_genes.txt',head=F,sep='\t')
VS16_target <- merge(VS16,GO_target,by.x='V1',by.y='V1')
WBC3_target <- merge(WBC3,GO_target,by.x='V1',by.y='V1')
AP13_target <- merge(AP13,GO_target,by.x='V1',by.y='V1')
DAC6_target <- merge(DAC6,GO_target,by.x='V1',by.y='V1')
setdiff(intersect(VS16_target[,1],DAC6_target[,1]),unique(AP13_target[,1],WBC3_target[,1]))
TPM_target <- TPM_target[setdiff(intersect(VS16_target[,1],DAC6_target[,1]),union(AP13_target[,1],WBC3_target[,1])),]
pdf('signal_transduction_genes_only_upland_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
TPM_target <- TPM[setdiff(intersect(VS16_target[,1],DAC6_target[,1]),unique(AP13_target[,1],WBC3_target[,1])),]
pdf('signal_transduction_genes_only_upland_TPM.pdf')
heatmap.2(as.matrix(TPM_target),trace="none",col = colorRampPalette(c("yellow","red"))(21),dendrogram='row',cex.axis=0.4,notecex=0.4,cexRow=0.4,cexCol=0.4,Rowv = TRUE,Colv=FALSE,lmat=rbind(4:3,2:1),lwid=c(0.8,4),lhei=c(0.8,4))
dev.off()
FC_target <- FC[setdiff(intersect(VS16_target[,1],DAC6_target[,1]),union(AP13_target[,1],WBC3_target[,1])),]

#######################################################################################
###############															###############
###############						PCA metabolites 					###############
###############															###############
#######################################################################################
setwd('D:\\Switchgrass_allelic_specific\\Collaborate_with_Anne-Sophie_Xinging\\Correlation_between_metabolome_and_transcriptome\\PCA_metabolites')
dat <- read.table('Metabolites_matrix_XX_exp1.txt', row.names=1,head=T,stringsAsFactors=F)
pca <- prcomp(t(dat), center = TRUE,scale. = TRUE)
summary(pca)
str(pca)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(pca)
ggbiplot(pca, labels=rownames(dat))
group <- c(rep("VS16_D1", 10), rep("VS16_D2",9), rep("WBC3_D1", 10),rep("WBC3_D2", 10), rep("F1_D1", 14), rep("F1_D2", 14))
pdf('Metabolites_matrix_XX_exp1_PCA.pdf')
ggbiplot(pca,ellipse=TRUE,  labels=NULL, groups=group, pc.biplot = FALSE, var.axes = FALSE)
dev.off()
pdf('Metabolites_matrix_XX_exp1_PCA_PC3_4.pdf')
ggbiplot(pca,ellipse=TRUE, choices=c(3,4), labels=NULL, groups=group, pc.biplot = FALSE, var.axes = FALSE)
dev.off()
pdf('Metabolites_matrix_XX_exp1_PCA_PC2_3.pdf')
ggbiplot(pca,ellipse=TRUE, choices=c(2,3), labels=NULL, groups=group, pc.biplot = FALSE, var.axes = FALSE)
dev.off()

# showing contributions of varibles to PCs
#install.packages('factoextra')
library(factoextra)
fviz_contrib(pca, choice = "var", axes = 1, top = 10)

# median values among replicates
library(matrixStats)
dat <- cbind(dat,'VS16_D1_median' = rowMedians(as.matrix(dat[,1:10])))
dat <- cbind(dat,'VS16_D2_median' = rowMedians(as.matrix(dat[,11:19])))
dat <- cbind(dat,'WBC3_D1_median' = rowMedians(as.matrix(dat[,20:29])))
dat <- cbind(dat,'WBC3_D2_median' = rowMedians(as.matrix(dat[,30:39])))
dat <- cbind(dat,'F1_D1_median' = rowMedians(as.matrix(dat[,40:53])))
dat <- cbind(dat,'F1_D2_median' = rowMedians(as.matrix(dat[,54:67])))
dat2 <- dat[,68:73]
pca <- prcomp(t(dat2), center = TRUE,scale. = TRUE)
group <- as.factor(colnames(dat2))
ggbiplot(pca,ellipse=TRUE,  labels=NULL, groups=group, pc.biplot = FALSE, var.axes = FALSE)



#################################################################################################
################                                                                 ################
################                       Analysis for Thilanka's data              ################
################                                                                 ################
#################################################################################################
setwd('D:\\Switchgrass_allelic_specific\\04_presentation\\Data_from_Thilanka')
h30 <- read.table("GO_30_m_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h1 <- read.table("GO_1_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h3 <- read.table("GO_3_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h6 <- read.table("GO_6_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h16 <- read.table("GO_16_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h24 <- read.table("GO_24_U.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
all_go <- unique(c(h30[,1],h1[,1],h3[,1],h6[,1],h16[,1],h24[,1]))
res <- matrix(data = NA, nrow=length(all_go),ncol=6)
rownames(res) <- all_go
colnames(res) <- c('h30','h1','h3','h6','h16','h24')
for(GO in all_go){
	if(length(h30[h30[,1]==GO,2]) > 0) res[GO,'h30'] <- h30[h30[,1]==GO,2]
	if(length(h1[h1[,1]==GO,2]) > 0) res[GO,'h1'] <- h1[h1[,1]==GO,2]
	if(length(h3[h3[,1]==GO,2]) > 0) res[GO,'h3'] <- h3[h3[,1]==GO,2]
	if(length(h6[h6[,1]==GO,2]) > 0) res[GO,'h6'] <- h6[h6[,1]==GO,2]
	if(length(h16[h16[,1]==GO,2]) > 0) res[GO,'h16'] <- h16[h16[,1]==GO,2]
	if(length(h24[h24[,1]==GO,2]) > 0) res[GO,'h24'] <- h24[h24[,1]==GO,2]
	}

res <- -log10(res)
res[is.na(res)] <- 0
res <- res[nrow(res):1,]
res[res > 10] <- 10
library(pvclust)
library(gplots)
pdf('Enriched_GO_up-regulation_all_6_timepoints.pdf')
lmat <- rbind(4:3,2:1)
lwid = c(0.8,4)
lhei = c(0.8,4)
h <- heatmap.2(res,trace="none",col = colorRampPalette(c("white","red"))(8),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.1,Rowv = FALSE,Colv=FALSE,lmat=lmat,lwid=lwid,lhei=lhei)
dev.off()

unique(c(h30[,2],h1[,2],h3[,2],h6[,2],h16[,2],h24[,2]))



























