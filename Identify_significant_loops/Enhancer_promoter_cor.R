##### Calculate enhancers and genes coorelations of all Enhancer-Promoter loops
##### author: Zhen Y

##### usage: prepare config file contains all parameters that are needed
##### parse config parameters to R script: Enhancer_promoter_cor.R
##### config file is a yaml formated file, with example below:

# config.yml
default:
	expr: ./expression_matrix.csv
	acSig: ./enhancers_signal.csv
	group: Group3
	protein: ./protein_coding_genes.bed
	phenotype: ./phenotype.csv
	cutoff: 0.1
	annotation: ./Enhancer_promoter_loops.bed

#####  explanation:
#####  expression_matrix.csv： expression matrix of normalized genes counts of all samples; row: genes; column: samples
#####  enhancers_signal.csv： intensity matrix of normalized enhancer signal of all samples; row: enhancers; column: samples
#####  protein_coding_genes.bed: bed files with protein coding genes list
#####  phenotype.csv: phenotype file with row as samples, column as meta information
#####  cutoff: significant FDR cutoff
#####  Enhancer_promoter_loops.bed: enhancer-promoter pairs annotated by left and right anchors




########  R script Enhancer_promoter_cor.R

library(tidyverse)
library(qvalue)
options(stringsAsFactors=F)


# parse parameters from yml
pars = config::get(file = "config.yml")
expr = pars$expr
acSig = pars$acSig
group = pars$group
protein = pars$protein
pheno = pars$phenotype
cutoff = pars$cutoff
anno = pars$annotation

# loading files
cpm=read.csv(expr, check.names=F, row.names = "ID")
cpm = as.matrix(cpm)
ac=read.csv(acSig, check.names=F, row.names = "ID")
ac = as.matrix(ac)
pro = read.table(protein)
phe = read.csv(pheno)
anno = read.table(anno, sep="\t")



################################
################################ 	First part
################################

samples = phe$sample[which(phe$group == group)] ## subseting samples by group
## subset ac signal matrix by enhancer ID at anno file column 4
ac1 = ac[rownames(ac) %in% anno$V4, colnames(ac) %in% samples]; ac1 = as.matrix(ac1) 
rbe = cpm[, colnames(cpm) %in% samples]; rbe = as.matrix(rbe)

nrow(ac1) == nrow(anno)
Hugelist = list()
for (i in seq(nrow(anno))) {

	# enhancer-promoter pairs is based on enhancer ID
	# column 7 contain potential target genes of each enhancer that validated by loops, separated by ; 
	genes = unique(unlist(strsplit(anno$V7[i],";"))); genes = genes[genes %in% pro$V4]
	# make sure target genes are protein coding
	if (length(genes) > 0) {
	m1 = rbe[rownames(rbe) %in% genes, ]; m2 = ac1[rownames(ac1) %in% anno$V4[i],]
	
	if (nrow(m1) > 1) {
		res = apply(m1, 1, function(d) {
			res1 = cor.test(d, m2, method = "spearman") # spearman correlation 
			return(c(res1$p.value, res1$estimate))
			})

		res = t(res); 
		rownames(res) = genes; 
		colnames(res) = c("", "rho")
		res = res[order(res[,1]), ]
		Hugelist[[anno$V4[i]]] = res

	} else {

		res = cor.test(m1, m2, method = "spearman")
		res = c(res$p.value, res$estimate)
		res = t(matrix(res)); 
		rownames(res) = genes; 
		colnames(res) = c("", "rho")
		Hugelist[[anno$V4[i]]] = res
	}}
}

save(file = paste0(group,"_enhancer-promoter.Rdata"), Hugelist)


################################
################################ 	Second part	

MM = do.call(rbind, Hugelist); 
# assign enhancer ID
col1 = rep(names(Hugelist), sapply(Hugelist, nrow)); 
genes = rownames(MM)
MM = data.frame(MM); 
MM$genes = genes; 
MM$enhancer = col1; 
# only keep loops with positive coorelation
MM = MM[which(MM$rho>0), ]
# calculate FDR by package qvalue
colnames(MM)[1] = "pvalue"
FDR1 = qvalue(as.numeric(MM$pvalue))
MM$fdr = FDR1
# saving all loops
write.csv(MM, paste0(group, "_enhancers-promoters.txt"), sep="\t", quote = F, row.names = F)

MM = MM[which(as.numeric(MM$fdr) < cutoff), ]
write.csv(MM, paste0(group, "_", cutoff, "_enhancers-promoters.txt"), quote = F, sep = "\t", row.names = F)



