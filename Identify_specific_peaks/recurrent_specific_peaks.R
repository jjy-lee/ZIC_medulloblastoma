##### Recurrent specific peaks
##### author: Zhen Y

##### identify recurrent group-specific peaks


options(stringsAsFactors=F)
library(tidyverse)


# recurrence table for each group
# file contain two columns:
# first column: ID: peak and overlap element pair after annotation
# second column: Recur: recurrence of peak element pair in each group
g3 = read.csv("Group3_recurrence_matrix.csv", check.names = F)
g4 = read.csv("Group4_recurrence_matrix.csv", check.names = F)
shh = read.csv("SHH_recurrence_matrix.csv", check.names = F)
wnt = read.csv("WNT_recurrence_matrix.csv", check.names = F)

# merge four groups
x = Reduce(function(x, y) full_join(x = x, y = y, by = "ID"), list(wnt,shh,g3,g4)); x = as.data.frame(x)
names(x)[2:5] = c("WNT", "SHH", "G3", "G4")
rownames(x) = x$ID; x = x[,-1]
x = as.matrix(x); x[which(is.na(x))] = 0

# recurrent specific
# WNT group level 1
WNT = rownames(x)[which(x[,1] >= 3 & x[,2] < 2 & x[,3] < 2 & x[,4] < 2)]; 
# WNT group level 2
WNT1 = rownames(x)[which(x[,1] >= 3 & x[,3] < 2 & x[,4] < 2)]; 

# SHH group level 1
SHH = rownames(x)[which(x[,2] >= 3 & x[,1] < 2 & x[,3] < 2 & x[,4] < 2)]
# SHH group level 2
SHH1 = rownames(x)[which(x[,2] >= 3 & x[,3] < 2 & x[,4] < 2)]; 

# G3 group level 1
G3 = rownames(x)[which(x[,3] >= 3 & x[,1] < 2 & x[,2] < 2 & x[,4] < 2)];
# G3 group level 2
G31 = rownames(x)[which(x[,3] >= 3 & x[,1] < 2 & x[,2] < 2)]; 

# G4 group level 1
G4 = rownames(x)[which(x[,4] >= 3 & x[,1] < 2 & x[,2] < 2 & x[,3] < 2)]; 
# G4 group level 2
G41 = rownames(x)[which(x[,4] >= 3 & x[,1] < 2 & x[,2] < 2)]; 

WNTSHH = setdiff(WNT1, WNT)
SHHWNT = setdiff(SHH1, SHH)
G34 = setdiff(G31, G3)
G43 = setdiff(G41, G4)
# organize
df = data.frame(peak = c(WNT,SHH,G3,G4,WNTSHH, SHHWNT, G34, G43), 
	group = rep(c("WNT", "SHH", "Group3", "Group4", "WNT/SHH", "SHH/WNT", "Group3/4", "Group4/3"), 
		c(length(WNT), length(SHH), length(G3), length(G4), length(WNTSHH), length(SHHWNT), length(G34), length(G43))))

# assign recurrence to new data frame
x = as.data.frame(x)
df$WNT_recur = x$WNT[match(df$peak, rownames(x))]
df$SHH_recur = x$SHH[match(df$peak, rownames(x))]
df$G3_recur = x$G3[match(df$peak, rownames(x))]
df$G4_recur = x$G4[match(df$peak, rownames(x))]
df1 = df %>% group_by(peak) %>% mutate(group=paste(unique(group), collapse=","))
df1 = df1[!duplicated(df1$peak), ]; df = as.data.frame(df1)

write.table(df, "recurrent_specific_peaks.xls", sep="\t", quote=F, col.names=T, row.names=F)
