##### signal specific peaks
##### author: Zhen Y

##### identify signal intensity group-specific peaks



options(stringsAsFactors=F)
library(tidyverse)
# read significantly up and down regulated peaks from DESeq2 results
# file contain two column: peak ID and UP/DOWN regulated (e.g. UP in G3_G4 means this peak up-regulated in G3 compare to G4)

# Group3-specific 
##  UP type1
up1=read.table("G3_G4_UP.txt")$V1; up2=read.table("G3_SHH_UP.txt")$V1; up3=read.table("G3_WNT_UP.txt")$V1
g3_spe_up = intersect(up1, intersect(up2, up3));


## DOWN  type1
down1=read.table("G3_G4_DOWN.txt")$V1; down2=read.table("G3_SHH_DOWN.txt")$V1; down3=read.table("G3_WNT_DOWN.txt")$V1
g3_spe_down = intersect(down1, intersect(down2, down3))


write.table(g3_spe_up, "Group3_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(g3_spe_down, "Group3_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

## UP type2
up1=read.table("G3_SHH_UP.txt")$V1; up2=read.table("G3_WNT_UP.txt")$V1
g3_spe_up1 = intersect(up1, up2);
g3_spe_up1 = setdiff(g3_spe_up1, g3_spe_up)

## DOWN  type2
down1=read.table("G3_SHH_DOWN.txt")$V1; down2=read.table("G3_WNT_DOWN.txt")$V1
g3_spe_down1 = intersect(down1, down2)
g3_spe_down1 = setdiff(g3_spe_down1, g3_spe_down)

write.table(g3_spe_up1, "Group3_2_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(g3_spe_down1, "Group3_2_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

# Group4
## UP type1
up1=read.table("G3_G4_DOWN.txt")$V1; up2=read.table("G4_SHH_UP.txt")$V1; up3=read.table("G4_WNT_UP.txt")$V1
g4_spe_up = intersect(up1, intersect(up2, up3));


## DOWN  type1
down1=read.table("G3_G4_UP.txt")$V1; down2=read.table("G4_SHH_DOWN.txt")$V1; down3=read.table("G4_WNT_DOWN.txt")$V1
g4_spe_down = intersect(down1, intersect(down2, down3))


write.table(g4_spe_up, "Group4_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(g4_spe_down, "Group4_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

## UP type2
up1=read.table("G4_SHH_UP.txt")$V1; up2=read.table("G4_WNT_UP.txt")$V1
g4_spe_up1 = intersect(up1, up2);
g4_spe_up1 = setdiff(g4_spe_up1, g4_spe_up)

## DOWN  type2
down1=read.table("G4_SHH_DOWN.txt")$V1; down2=read.table("G4_WNT_DOWN.txt")$V1
g4_spe_down1 = intersect(down1, down2)
g4_spe_down1 = setdiff(g4_spe_down1, g4_spe_down)

write.table(g4_spe_up1, "Group4_2_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(g4_spe_down1, "Group4_2_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

#SHH
## UP type1
up1=read.table("SHH_WNT_UP.txt")$V1; up2=read.table("G4_SHH_DOWN.txt")$V1; up3=read.table("G3_SHH_DOWN.txt")$V1
sh_spe_up = intersect(up1, intersect(up2, up3));


## DOWN  type1
down1=read.table("SHH_WNT_DOWN.txt")$V1; down2=read.table("G4_SHH_UP.txt")$V1; down3=read.table("G3_SHH_UP.txt")$V1
sh_spe_down = intersect(down1, intersect(down2, down3))


write.table(sh_spe_up, "SHH_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(sh_spe_down, "SHH_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

## UP type2
up1=read.table("G4_SHH_DOWN.txt")$V1; up2=read.table("G3_SHH_DOWN.txt")$V1
sh_spe_up1 = intersect(up1, up2);
sh_spe_up1 = setdiff(sh_spe_up1, sh_spe_up)

## DOWN  type2
down1=read.table("G4_SHH_UP.txt")$V1; down2=read.table("G3_SHH_UP.txt")$V1
sh_spe_down1 = intersect(down1, down2)
sh_spe_down1 = setdiff(sh_spe_down1, sh_spe_down)

write.table(sh_spe_up1, "SHH_2_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(sh_spe_down1, "SHH_2_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

# WNT
## UP type1
up1=read.table("SHH_WNT_DOWN.txt")$V1; up2=read.table("G4_WNT_DOWN.txt")$V1; up3=read.table("G3_WNT_DOWN.txt")$V1
wn_spe_up = intersect(up1, intersect(up2, up3));


## DOWN  type1
down1=read.table("SHH_WNT_UP.txt")$V1; down2=read.table("G4_WNT_UP.txt")$V1; down3=read.table("G3_WNT_UP.txt")$V1
wn_spe_down = intersect(down1, intersect(down2, down3))


write.table(wn_spe_up, "WNT_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(wn_spe_down, "WNT_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)

## UP type2
up1=read.table("G4_WNT_DOWN.txt")$V1; up2=read.table("G3_WNT_DOWN.txt")$V1
wn_spe_up1 = intersect(up1, up2);
wn_spe_up1 = setdiff(wn_spe_up1, wn_spe_up)

## DOWN  type2
down1=read.table("G4_WNT_UP.txt")$V1; down2=read.table("G3_WNT_UP.txt")$V1
wn_spe_down1 = intersect(down1, down2)
wn_spe_down1 = setdiff(wn_spe_down1, wn_spe_down)

write.table(wn_spe_up1, "WNT_2_specific_peaks_UP.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(wn_spe_down1, "WNT_2_specific_peaks_DOWN.txt", sep="\t", quote=F, col.names=F, row.names=F)


dfup = data.frame(peak = c(wn_spe_up,sh_spe_up,g3_spe_up,g4_spe_up,wn_spe_up1,sh_spe_up1,g3_spe_up1,g4_spe_up1), 
	group = rep(c("WNT", "SHH", "Group3", "Group4", "WNT/SHH", "SHH/WNT", "Group3/4", "Group4/3"), 
		c(length(wn_spe_up), length(sh_spe_up), length(g3_spe_up), length(g4_spe_up), length(wn_spe_up1), 
			length(sh_spe_up1), length(g3_spe_up1), length(g4_spe_up1))))

dfdown = data.frame(peak = c(wn_spe_down,sh_spe_down,g3_spe_down,g4_spe_down,wn_spe_down1,sh_spe_down1,g3_spe_down1,g4_spe_down1), 
	group = rep(c("WNT", "SHH", "Group3", "Group4", "WNT/SHH", "SHH/WNT", "Group3/4", "Group4/3"), 
		c(length(wn_spe_down), length(sh_spe_down), length(g3_spe_down), length(g4_spe_down), length(wn_spe_down1), 
			length(sh_spe_down1), length(g3_spe_down1), length(g4_spe_down1))))


write.table(dfup, "specific_upreg_peaks.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(dfdown, "specific_downreg_peaks.txt", sep="\t", quote=F, col.names=F, row.names=F)


