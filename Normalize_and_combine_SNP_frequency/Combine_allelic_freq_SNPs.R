#Zhen You
####################
####################
####################		Combine allelic frequencies of sequential het. SNPs inside a functional element
####################
####################
#
#  Project: ZIC1 is a context dependent cancer driver in the rhombic lip
#  allelic frequencies of each SNP in each sample was calculated and stored in independent file
#  sample information sheet is needed
#  phased SNPs for each sample are needed
#  This script will output the aggregated allelic frequencies for genes/peaks for each sample based on phased SNPs allele ratios,
#  which are weighted by the reads depth

#  usage: Combine_allelic_freq_SNPs.R subgroup cutoff output_dir
#  subgroup: WNT, SHH, Group3, Group4
#  cutoff: read depth cutoff
#  output_dir: folder to store result table

#  SNP information needs to be prepared in _results.txt;
#  _results.txt contain 10 columns: chr, pos, ref, alt, total depth, alt depth, imbalance ratio, p-value, peaks/genes, phase info

options(stringsAsFactors=F)
library(tidyverse)
library(data.table)
library(stringr)
library(metap)

args = commandArgs(T)
subgroup = args[1]
cutoff = args[2]
output = args[3]

saminfo = read.csv("SampleInfo_metasheet.csv"); ## sample information
saminfo = saminfo$sample[which(saminfo$group==subgroup)]; files = paste0(saminfo, "_results.txt");  # allelic frequencies results of each sample

for (q in seq(length(files))) {
	print(q)
	x=read.table(files[q]); x = x[which(x$V5>=cutoff), ];
	## total depth cutoff set to 5 for ChIP-Seq, 10 for RNA-Seq
	x = x[!duplicated(x[,c(1:6,9)]),]
	# remove redundant SNPs
	x$V11 = p.adjust(x$V8, method="fdr");
	x$ID = paste(x$V1, x$V2, x$V9, sep = "_");
	# make unique ID chr_pos_peaks/genes

	x = x %>% group_by(V9) %>% mutate(V5 = V5/sum(V5))
	# group by peak ID or gene
	# generate weights for each SNP in a element

	x = as.data.frame(x); targets = unique(x$V9) ## targets contain peaks/genes
	element = c(); ratios = c(); mfdr = c();
	for (t in targets) {
		mdf = x[which(x$V9 == t), ]
		# haplotype of each SNP extract from phasing result, e.g., 1|0 or 0|1
		# haplotype information was convert to haplotype 0 and haplotype 1 in _results.txt
		pp = sumz(mdf$V11, weights = mdf$V5); pp = as.numeric(pp$p)
		# multiple fdr were integrated using metap package sumz function, and weighted by total depth

		hp0 = ifelse(mdf$V10 == 0, mdf$V7, 1 - mdf$V7)
		hp1 = 1 - hp0
		### calculate allele ratios on each haplotype
		hp_a1 = sum((hp1 - 0.5) * mdf$V5); hp_a0 = sum((hp0 - 0.5) * mdf$V5);
		## transform allele frequencies (positive ratio for dominant haplotype) 
		## allele ratios from both haplotypes were weighted by read depth
		hp_a = abs(max(hp_a1, hp_a0))
		
		element = c(element, t); ratios = c(ratios, hp_a); mfdr = c(pp)
	}
	wdf = data.frame(element = element, ratios = ratios, fdr = mfdr)
	ID = gsub("_results.txt", "", basename(files[q])
	write.table(wdf, paste0(output, "/", ID, "_output.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
	# write down tables with peaks/genes, allele frequencies, FDR
}
