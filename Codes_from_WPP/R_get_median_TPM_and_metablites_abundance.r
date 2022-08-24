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
