##### CRC analysis: Enrichment score of TFs calculated from CRC program results 
##### author: Zhen Y

##### This script calculate TF enrichment score from auto-regulatory circuits generated from CRC program
##### This script also calculate group-specific master TFs by their enrichment scores
##### Usage: Rscript CRC_enrichment_score.R path_to_CRC_results metadata.csv
##### path_to_CRC_results: path to where CRC results stored
##### metadata.csv: sample information, must contain two column: group: subgroups, sample: sample ID

options(stringsAsFactors=F)
library(dplyr)
args = commandArgs(T)

if (length(args) != 3) {
	cat("Error: arguments incorrect;\nUsage: Rscript CRC_enrichment_score.R path_to_CRC_results metadata.csv\n")
	quit("no", 1, FALSE)
}

dirs = args[1]
metas = args[2]

info=read.csv(metas, colClass = "character")
reslist = vector("list", 4)
names(reslist) = c("WNT","SHH","Group3", "Group4")

for (group in c("WNT", "SHH", "Group3", "Group4")) {
	samples = info$sample[which(info$group == group)]
	# file structure:
	# assume CRC results of each sample are stored at path_to_CRC_results, with each sample has its own folder, named as sample ID
	cliques = paste0(dirs, "/", samples,"/",samples,"_EDGE_LIST.txt")
	# *_EDGE_LIST.txt: List of edges describing which TFs interact with each other.
	auto = paste0(dirs, "/", samples,"/",samples,"_SELF_LOOPS.txt")
	# *_SELF_LOOPS.txt: auto-regulatory TFs list
	# file name: use sample ID as analysis name when running CRC program

	# reading files
	cli = lapply(cliques, function(i) {
		x1 = read.table(i, sep = "\t", header = T)
	})

	aut = lapply(auto, function(i) {
		x1 = read.table(i, sep = "\t")$V1;
	})

	for (i in seq(length(cli))) {
		dff = data.frame(); 

		for (l in aut[[i]]) {
			a1 = intersect(cli[[i]]$To[which(cli[[i]]$From==l)], cli[[i]]$From[which(cli[[i]]$To==l)]) ## make sure they are loops (both target at and targeted by)
			a2 = intersect(aut[[i]], a1)  ## TFs in loops must be self-loops
			df = data.frame(V1 = rep(l, length(a2)), V2 = a2); 
			dff = rbind(dff, df)
		}

		write.table(dff, gsub("EDGE_LIST", "auto-regulatory", cliques[i]), sep = "\t", quote = F, col.names = F, row.names = F)
	}
	# save auto-regulatory TF pairs for each sample

	#### calculate the percentage matrix for each sample CRC score rank matrix
	cliques = paste(dirs, "/", samples,"/",samples,"_auto-regulatory.txt",sep="")

	cli = lapply(cliques, function(i) {
		x1 = read.table(i, sep = "\t")
	})

	mylist = vector("list", length(cli))

	for (i in seq(length(cli))) {

		candidates = lapply(as.character(unique(cli[[i]]$V1)), function(l) {
			return(length(unique(cli[[i]]$V2[grep(l, cli[[i]]$V1)])))
		})
		# For each auto-regulatory TFj, This step return number of TFs that target at/targeted by TFj (inter-connect with TFj)
		perc = sapply(candidates, function(p) length(p)/length(unique(cli[[i]]$V2)))
		# calculate percentage
		names(perc) = as.character(unique(cli[[i]]$V1))
		mylist[[i]] = perc
	}

	names(mylist) = samples
	### sort percentage matrix 
	x = Reduce(function(x,y) merge(x = x, y = y, by = "V1", all=T), lapply(mylist, function(i) {p=data.frame(i); p$V1=rownames(p); return(p)}))
	
	names(x) = c("TF",samples); x = x[which(rowSums(is.na(x)==F)>2), ]  # TF recur at least in two samples
	reslist[[group]] = x
}


fallp = Reduce(function(x,y) merge(x = x, y = y, by = "TF", all=T), reslist)
colnames(fallp) = c("TF", paste(rep(c("WNT","SHH","G3","G4"), sapply(reslist, length)-1), colnames(fallp)[2:ncol(fallp)], sep = "_"))
write.csv(fallp, "enrichment_score_perc_matrix.csv", quote = F, row.names=F)


# calculate group-specific TFs
rownames(fallp) = fallp$TF; fallp = fallp[,-1]
fallp = as.matrix(fallp); 
fallp[which(is.na(fallp))] = 0
Matrix = matrix(0, nrow(fallp), 6)

for (i in seq(nrow(fallp))) {
	p = stack(fallp[i,]); p$ind = gsub("\\.\\d*", "", p$ind) # create data frame with one column as groups and the other column as percentages

	p$ind = factor(p$ind); 
	# calculating p value from t-test
	Matrix[i,1] = t.test(p$values[which(p$ind=="G3")], p$values[which(p$ind=="WNT")])$p.value
	Matrix[i,2] = t.test(p$values[which(p$ind=="G3")], p$values[which(p$ind=="SHH")])$p.value
	Matrix[i,3] = t.test(p$values[which(p$ind=="G3")], p$values[which(p$ind=="G4")])$p.value
	Matrix[i,4] = t.test(p$values[which(p$ind=="WNT")], p$values[which(p$ind=="G4")])$p.value
	Matrix[i,5] = t.test(p$values[which(p$ind=="SHH")], p$values[which(p$ind=="G4")])$p.value
	Matrix[i,6] = t.test(p$values[which(p$ind=="WNT")], p$values[which(p$ind=="SHH")])$p.value	
}

# matrix contains p values of all six comparisons
colnames(Matrix) = c("G3-WNT", "G3-SHH", "G3-G4", "G4-WNT", "G4-SHH", "WNT-SHH")	
rownames(Matrix) = rownames(fallp)

# 1, 4, 6 WNT 
# 2, 5, 6 SHH 
# 1, 2, 3 G3 
# 3, 4, 5 G4 
WNT_spe = c()
SHH_spe = c()
G3_spe = c()
G4_spe = c()

for (i in seq(nrow(Matrix))) {
	if (all(c(1,4,6) %in% which(Matrix[i,]<0.05))) {  # p value < 0.05 compare to three other groups
		WNT_spe = c(WNT_spe, rownames(Matrix)[i])
	} 
	if (all(c(2,5,6) %in% which(Matrix[i,]<0.05))) {
		SHH_spe = c(SHH_spe, rownames(Matrix)[i])
	} 
	if (all(c(1,2,3) %in% which(Matrix[i,]<0.05))) {
		G3_spe = c(G3_spe, rownames(Matrix)[i])
	} 
	if (all(c(3,4,5) %in% which(Matrix[i,]<0.05))) {
		G4_spe = c(G4_spe, rownames(Matrix)[i])
	}
}

# G3/4 specifc 1, 2, 4, 5 < 0.05; 3 > 0.05
# WNT/SHH specifc 1, 2, 4, 5 < 0.05; 6 > 0.05
G34_spe = c()
wntshh_spe = c()

for (i in seq(nrow(Matrix))) {
	if (all(c(1,2,4,5) %in% which(Matrix[i,]<0.05)) & Matrix[i,6]>0.05) {
		wntshh_spe = c(wntshh_spe, rownames(Matrix)[i])
	} 
	if (all(c(1,2,4,5) %in% which(Matrix[i,]<0.05)) & Matrix[i,3]>0.05) {
		G34_spe = c(G34_spe, rownames(Matrix)[i])
	}
}

### Group specific auto-regulatory TFs
specific_auto = data.frame(c(WNT_spe, SHH_spe, G3_spe, G4_spe, wntshh_spe, G34_spe), 
	rep(c("WNT","SHH","G3","G4","WNT/SHH","G3/4"), sapply(list(WNT_spe, SHH_spe, G3_spe, G4_spe, wntshh_spe, G34_spe), length)))
names(specific_auto) = c("TF", "group")
specific_auto = specific_auto %>% group_by(TF) %>% mutate(group=paste(group, collapse=","))
write.table(specific_auto, "enrichment_score_group_specific.txt", sep="\t", quote=F, col.names=T, row.names=F)

# ranking group-specific TFs by average enrichment score
x = fallp[rownames(fallp) %in% specific_auto$TF, ]
WNT = apply(x[,grep("WNT", colnames(fallp))], 1, mean, na.rm=T)
SHH = apply(x[,grep("SHH", colnames(fallp))], 1, mean, na.rm=T)
G3 = apply(x[,grep("G3", colnames(fallp))], 1, mean, na.rm=T)
G4 = apply(x[,grep("G4", colnames(fallp))], 1, mean, na.rm=T)
df = data.frame(WNT = WNT, SHH = SHH, G3 = G3, G4 = G4); rownames(df) = rownames(x)
df$group = specific_auto$group[match(rownames(df), specific_auto$TF)]
write.table(df, "enrichment_score_group_specific_ranking.txt", sep="\t", quote=F, col.names=T, row.names=T)
