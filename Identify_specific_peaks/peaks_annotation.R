
##### Peaks annotation
##### author: Zhen Y

##### usage: Rscript peaks_annotation.R input.bed output.bed
##### input.bed: H3K27Ac or H3K27me3 peaks bed file (5 column bed format)
##### output.bed: output annotation file name
##### Create peaks annotations for both merged and individual peak file
##### Create peaks genes element pairs for recurrence specific peaks identification 


options(stringsAsFactors=F)

args = commandArgs(T)

if (length(args) != 2) {
	cat("Error: arguments incorrect;\nUsage: Rscript peaks_annotation.R input.bed output.bed\n")
	quit("no", 1, FALSE)
}

# annotation of each bed file 
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(tools)
library(dplyr)

input = args[1]
output = args[2]

# Just regular five column bed format
peaks = read.table(input)

names(peaks)=c("chr","start","end","id","value","strand")
gr1 <- toGRanges(peaks, format="BED", header=FALSE)
ucsc.hg19.knownGene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
macs.anno <- annotatePeakInBatch(gr1, 
                                 AnnotationData=ucsc.hg19.knownGene, output=c("both"), FeatureLocForDistance="TSS")

macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                        orgAnn="org.Hs.eg.db", 
                        feature_id_type="entrez_id",
                        IDs2Add="symbol")
macs = data.frame(macs.anno)

write.table(macs, output, sep = "\t", quote=F, col.names = T, row.names = F)

# neat version of annotation
# file name suffix
suffix1 = file_ext(output)
# file name prefix
prefix1 = file_path_sans_ext(output)
prefix1 = paste(prefix1, "neatAnno", sep="_")
newout = paste(prefix1, suffix1, sep = ".")
macs_neat = macs[,c(1:3,6,13,17)]
# multiple genes that overlap with enhancer were grouped and separated by ","
macs_neat = macs_neat %>% group_by(id) %>% mutate_at(c(5,6), paste0, collapse = ",")
macs_neat = data.frame(macs_neat); macs_neat = macs_neat[!duplicated(macs_neat$id), ]
write.table(macs_neat, macs_neat, sep = "\t", quote=F, col.names = T, row.names = F)



# classification of annotation for recurrent group-specific peaks
# priority promoter > exon > intron
# means if a peak overlaps with promoter and 5' UTR or first exon, this peak will be noted as promoter peak
# if a peak overlaps with a exon and the intron flanking that exon, this peak will be noted as exon peak

ucsc.hg19.exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)


macs$newanno = ""
candis = c("includeFeature", "overlapStart") 
# overlapStart: peaks overlap with TSS
# includeFeature: tht gene are locate inside the peak
# peak with these classification are regard as promoter peaks
macs$newanno[which(macs$insideFeature %in% candis)] = "promoter"
# Also include peaks with distance shorter than 2000 to TSS
# make sure the peak is not at the 3 prime of the gene
macs$newanno[which(abs(macs$distancetoFeature)<2000 & macs$insideFeature != "downstream")] = "promoter"
# distancetoFeature: distance to TSS

# intergenic regions are then selected out by peaks that do not overlap with genes and promoter regions
macs$newanno[which(macs$newanno != "promoter" & macs$insideFeature %in% c("downstream","upstream"))] = "intergenic"
macs$id = paste(macs$id, macs$feature, sep="_"); 

prom = macs[which(macs$newanno == "promoter"), ]
inter = macs[which(macs$newanno == "intergenic"), ]
write.table(prom, "promoter_peaks_Anno.bed", sep="\t", quote=F, col.names=T, row.names=F)
write.table(inter, "intergenic_peaks_Anno.bed", sep="\t", quote=F, col.names=T, row.names=F)

macs1 = macs[which(macs$newanno == ""), ];

# macs1 are remaining peaks overlap with gene body
gr1 <- toGRanges(macs1, format="BED", header=FALSE)
# find peaks overlap with exons
macs.anno <- annotatePeakInBatch(gr1, AnnotationData=ucsc.hg19.exons, output=c("both"), FeatureLocForDistance="TSS")
macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                    orgAnn="org.Hs.eg.db", 
                    feature_id_type="entrez_id",
                    IDs2Add="symbol")
macs = data.frame(macs.anno)

# peak overlaps with exon (both 5 and 3 prime), peak locate within a exon, and peak includes both exons and introns are regard as exon peak
macs$newanno = ""; macs$newanno[which(macs$insideFeature %in% 
	c("includeFeature","inside","overlapEnd","overlapStart"))] = "exon"
macs$newanno[which(macs$newanno != "exon")] = "intron";

exon1 = macs[which(macs$newanno == "exon"), ]
intron1 = macs[which(macs$newanno == "intron"), ]
write.table(exon1, "exon_peaks_Anno.bed", sep="\t", quote=F, col.names=T, row.names=F)
write.table(intron1, "intron_peaks_Anno.bed", sep="\t", quote=F, col.names=T, row.names=F)


