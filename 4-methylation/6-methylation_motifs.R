##### This script probably isn't the easiest way of doing this, but the #####
##### script produces a file that contains all of the different methylation #####
##### contexts separated. ####
##### This file can then be used as input in smoothing #####

###Loading all clones###
All_clones <- read.table("/path/to/filtered/covs/All_chromosomes.bismark.cov", header=T)
names(All_clones) <- c("V1", "V2", "swasp001_meth", "swasp005_meth", "swasp044_meth", "swasp046_meth", "swasp113_meth", "swasp116_meth", "mean", "std_dev")



###Loading motifs###
# Load in swasp001
CHG_swaspsaevar1_1_14 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar1-1-14_S88_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar1_1_14 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar1-1-14_S88_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar1_1_14 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar1-1-14_S88_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp001<-merge(CpG_swaspsaevar1_1_14,All_clones[,c("V1","V2","swasp001_meth")],by=c("V1","V2"))
CHG_swasp001<-merge(CHG_swaspsaevar1_1_14,All_clones[,c("V1","V2","swasp001_meth")],by=c("V1","V2"))
CHH_swasp001<-merge(CHH_swaspsaevar1_1_14,All_clones[,c("V1","V2","swasp001_meth")],by=c("V1","V2"))

swasp001<-merge(CpG_swasp001,CHG_swasp001,by=c("V1","V2"), all=T)
swasp001<-merge(swasp001,CHH_swasp001,by=c("V1","V2"), all=T)

names(swasp001) <- c("chromosome","position","CpG_swasp001","CHG_swasp001","CHH_swasp001")

rm(CHG_swaspsaevar1_1_14,CHH_swaspsaevar1_1_14,CpG_swaspsaevar1_1_14,CpG_swasp001,CHG_swasp001,CHH_swasp001)


# Load in swasp005
CHG_swaspsaevar5_8_17 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar5-8-17_S91_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar5_8_17 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar5-8-17_S91_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar5_8_17 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar5-8-17_S91_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp005<-merge(CpG_swaspsaevar5_8_17,All_clones[,c("V1","V2","swasp005_meth")],by=c("V1","V2"))
CHG_swasp005<-merge(CHG_swaspsaevar5_8_17,All_clones[,c("V1","V2","swasp005_meth")],by=c("V1","V2"))
CHH_swasp005<-merge(CHH_swaspsaevar5_8_17,All_clones[,c("V1","V2","swasp005_meth")],by=c("V1","V2"))

swasp005<-merge(CpG_swasp005,CHG_swasp005,by=c("V1","V2"), all=T)
swasp005<-merge(swasp005,CHH_swasp005,by=c("V1","V2"), all=T)

names(swasp005) <- c("chromosome","position","CpG_swasp005","CHG_swasp005","CHH_swasp005")

rm(CHG_swaspsaevar5_8_17,CHH_swaspsaevar5_8_17,CpG_swaspsaevar5_8_17,CpG_swasp005,CHG_swasp005,CHH_swasp005)


# Load in swasp044
CHG_swaspsaevar44_3_22 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar44-3-22_S92_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar44_3_22 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar44-3-22_S92_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar44_3_22 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar44-3-22_S92_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp044<-merge(CpG_swaspsaevar44_3_22,All_clones[,c("V1","V2","swasp044_meth")],by=c("V1","V2"))
CHG_swasp044<-merge(CHG_swaspsaevar44_3_22,All_clones[,c("V1","V2","swasp044_meth")],by=c("V1","V2"))
CHH_swasp044<-merge(CHH_swaspsaevar44_3_22,All_clones[,c("V1","V2","swasp044_meth")],by=c("V1","V2"))

swasp044<-merge(CpG_swasp044,CHG_swasp044,by=c("V1","V2"), all=T)
swasp044<-merge(swasp044,CHH_swasp044,by=c("V1","V2"), all=T)

names(swasp044) <- c("chromosome","position","CpG_swasp044","CHG_swasp044","CHH_swasp044")

rm(CHG_swaspsaevar44_3_22,CHH_swaspsaevar44_3_22,CpG_swaspsaevar44_3_22,CpG_swasp044,CHG_swasp044,CHH_swasp044)

# Load in swasp046
CHG_swaspsaevar46_1_1 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar46-1-1_S94_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar46_1_1 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar46-1-1_S94_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar46_1_1 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar46-1-1_S94_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp046<-merge(CpG_swaspsaevar46_1_1,All_clones[,c("V1","V2","swasp046_meth")],by=c("V1","V2"))
CHG_swasp046<-merge(CHG_swaspsaevar46_1_1,All_clones[,c("V1","V2","swasp046_meth")],by=c("V1","V2"))
CHH_swasp046<-merge(CHH_swaspsaevar46_1_1,All_clones[,c("V1","V2","swasp046_meth")],by=c("V1","V2"))

swasp046<-merge(CpG_swasp046,CHG_swasp046,by=c("V1","V2"), all=T)
swasp046<-merge(swasp046,CHH_swasp046,by=c("V1","V2"), all=T)

names(swasp046) <- c("chromosome","position","CpG_swasp046","CHG_swasp046","CHH_swasp046")

rm(CHG_swaspsaevar46_1_1,CHH_swaspsaevar46_1_1,CpG_swaspsaevar46_1_1,CpG_swasp046,CHG_swasp046,CHH_swasp046)

# Load in swasp113
CHG_swaspsaevar113_9_2 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar113-9-2_S97_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar113_9_2 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar113-9-2_S97_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar113_9_2 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar113-9-2_S97_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp113<-merge(CpG_swaspsaevar113_9_2,All_clones[,c("V1","V2","swasp113_meth")],by=c("V1","V2"))
CHG_swasp113<-merge(CHG_swaspsaevar113_9_2,All_clones[,c("V1","V2","swasp113_meth")],by=c("V1","V2"))
CHH_swasp113<-merge(CHH_swaspsaevar113_9_2,All_clones[,c("V1","V2","swasp113_meth")],by=c("V1","V2"))

swasp113<-merge(CpG_swasp113,CHG_swasp113,by=c("V1","V2"), all=T)
swasp113<-merge(swasp113,CHH_swasp113,by=c("V1","V2"), all=T)

names(swasp113) <- c("chromosome","position","CpG_swasp113","CHG_swasp113","CHH_swasp113")

rm(CHG_swaspsaevar113_9_2,CHH_swaspsaevar113_9_2,CpG_swaspsaevar113_9_2,CpG_swasp113,CHG_swasp113,CHH_swasp113)

# Load in swasp116
CHG_swaspsaevar116_6_28 <- read.table("/path/to/motif/lists/CHG_context_PE.swaspsaevar116-6-28_S98_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CHH_swaspsaevar116_6_28 <- read.table("/path/to/motif/lists/CHH_context_PE.swaspsaevar116-6-28_S98_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)
CpG_swaspsaevar116_6_28 <- read.table("/path/to/motif/lists/CpG_context_PE.swaspsaevar116-6-28_S98_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.methylated_list.sorted.positions.chr.txt", header=F)

CpG_swasp116<-merge(CpG_swaspsaevar116_6_28,All_clones[,c("V1","V2","swasp116_meth")],by=c("V1","V2"))
CHG_swasp116<-merge(CHG_swaspsaevar116_6_28,All_clones[,c("V1","V2","swasp116_meth")],by=c("V1","V2"))
CHH_swasp116<-merge(CHH_swaspsaevar116_6_28,All_clones[,c("V1","V2","swasp116_meth")],by=c("V1","V2"))

swasp116<-merge(CpG_swasp116,CHG_swasp116,by=c("V1","V2"), all=T)
swasp116<-merge(swasp116,CHH_swasp116,by=c("V1","V2"), all=T)

names(swasp116) <- c("chromosome","position","CpG_swasp116","CHG_swasp116","CHH_swasp116")

rm(CHG_swaspsaevar116_6_28,CHH_swaspsaevar116_6_28,CpG_swaspsaevar116_6_28,CpG_swasp116,CHG_swasp116,CHH_swasp116)



#Merge files into one large file
all_clones_all_motifs<-merge(swasp001, swasp005, by=c("chromosome", "position"), all=T)
all_clones_all_motifs<-merge(all_clones_all_motifs, swasp044, by=c("chromosome", "position"), all=T)
all_clones_all_motifs<-merge(all_clones_all_motifs, swasp046, by=c("chromosome", "position"), all=T)
all_clones_all_motifs<-merge(all_clones_all_motifs, swasp113, by=c("chromosome", "position"), all=T)
all_clones_all_motifs<-merge(all_clones_all_motifs, swasp116, by=c("chromosome", "position"), all=T)

rm(All_clones, swasp001, swasp005, swasp044, swasp046, swasp113, swasp116)

write.table(all_clones_all_motifs, "/wherever/you/want/to/save/this/all_clones_all_motifs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image(file="/path/to/all_clones_all_motifs.Rdata")

### Smooth motifs ###
### This part is for the smoothed motifs, after running 7-methylation_motifs_windows.py and smoothing ###
# Calculates means and standard deviations for the smoothed windows #
all_clones_all_motifs <- read.csv("/path/to/windowed_motif/all_clones_all_motifs.windowed.1Mbp_window.txt", header=T, sep="\t")
all_clones_all_motifs$CpG_mean<-apply(all_clones_all_motifs[,c(4,7,10,13,16,19)],1,mean, na.rm=T)
all_clones_all_motifs$CpG_sd<-apply(all_clones_all_motifs[,c(4,7,10,13,16,19)],1,sd, na.rm=T)
all_clones_all_motifs$CHG_mean<-apply(all_clones_all_motifs[,c(5,8,11,14,17,20)],1,mean, na.rm=T)
all_clones_all_motifs$CHG_sd<-apply(all_clones_all_motifs[,c(5,8,11,14,17,20)],1,sd, na.rm=T)
all_clones_all_motifs$CHH_mean<-apply(all_clones_all_motifs[,c(6,9,12,15,18,21)],1,mean, na.rm=T)
all_clones_all_motifs$CHH_sd<-apply(all_clones_all_motifs[,c(6,9,12,15,18,21)],1,sd, na.rm=T)

write.table(all_clones_all_motifs, "/path/to/windowed_motif/all_clones_all_motifs.windowed.1Mbp_window.mean.txt", sep = "\t", quote = FALSE, row.names = FALSE)
