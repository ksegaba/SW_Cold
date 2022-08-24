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
