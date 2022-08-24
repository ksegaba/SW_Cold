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

