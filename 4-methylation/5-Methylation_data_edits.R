##### Script for making a file with all the filtered coverages #####


setwd("/path/to/filtered/covs/")

#Read SwAsp001
read.table("/path/to/filtered/covs/PE.swaspsaevar1-1-14_S88_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp001
names(SwAsp001)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")

#Read SwAsp005
read.table("/path/to/filtered/covs/PE.swaspsaevar5-8-17_S91_L006_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp005
names(SwAsp005)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")

#Read SwAsp113
read.table("/path/to/filtered/covs/PE.swaspsaevar113-9-2_S97_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp113
names(SwAsp113)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")

#Read SwAsp116
read.table("/path/to/filtered/covs/PE.swaspsaevar116-6-28_S98_L008_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp116
names(SwAsp116)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")

#Read SwAsp044
read.table("/path/to/filtered/covs/PE.swaspsaevar44-3-22_S92_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp044
names(SwAsp044)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")

#Read SwAsp046
read.table("/path/to/filtered/covs/PE.swaspsaevar46-1-1_S94_L007_R1_001_val_1_bismark_bt2_pe.deduplicated.chr.sorted.5-45_filter_cov.bismark.cov",head=F)->SwAsp046
names(SwAsp046)<-c("Scaffold","Start","End","%_methylation","Count_methylated","Count_unmethylated")


#Merge files into one large file
All_clones<-merge(SwAsp001[,c(1:2,4)],SwAsp005[,c(1:2,4)],by=c("Scaffold","Start"),all=T)
All_clones<-merge(All_clones,SwAsp044[,c(1:2,4)],by=c("Scaffold","Start"),all=T)
All_clones<-merge(All_clones,SwAsp046[,c(1:2,4)],by=c("Scaffold","Start"),all=T)
All_clones<-merge(All_clones,SwAsp113[,c(1:2,4)],by=c("Scaffold","Start"),all=T)
All_clones<-merge(All_clones,SwAsp116[,c(1:2,4)],by=c("Scaffold","Start"),all=T)
names(All_clones)<-c("Chromosome","Position","SwAsp001_meth%", "SwAsp005_meth%",
                     "SwAsp044_meth%", "SwAsp046_meth%", "SwAsp113_meth%", "SwAsp116_meth%")

#Remove the now unnecessary individual files
#rm(SwAsp001, SwAsp005, SwAsp044, SwAsp046, SwAsp113, SwAsp116)

#Calculate mean and standard deviation per position

All_clones$Mean<-apply(All_clones[,3:8],1,mean, na.rm=T)
All_clones$STD_dev<-apply(All_clones[,3:8],1,sd, na.rm=T)

#Write out into a table
write.table(All_clones, "/path/to/filtered/covs/All_chromosomes.bismark.cov", sep = "\t", quote = FALSE, row.names = FALSE)

#Release memory
rm(SwAsp001, SwAsp005, SwAsp044, SwAsp046, SwAsp113, SwAsp116)


### Just for checking everything is fine ###
#Calculate windowed averages
#window_size=50000
#na=data.frame()
#All_clones_windowed=data.frame()
#for (chromosome in unique(All_clones$Chromosome)) {
#  All_cloneschr <- All_clones[All_clones$Chromosome==chromosome,]
#  win_start = 0
#  win_end = 0
#  run=TRUE
#  while (run==TRUE) {
#    win_end=win_end+window_size
#    if (win_end < min(All_cloneschr$Position)) {
#      c=NA
#      na = append(na, c(chromosome, win_start, win_end, c,c,c,c,c,c))
#      win_start=win_end
#    }
#    if (length(All_cloneschr[win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end,])==0) {
#      c=NA
#      na = append(na, c(chromosome, win_start, win_end, c,c,c,c,c,c))
#      win_start=win_end
#    }
#    if (win_end > max(All_cloneschr$Position)) {
#      win_end = max(All_cloneschr$Position)
#      All_clones_windowed = rbind(All_clones_windowed, data.frame(chromosome, win_start, win_end,
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp001_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp005_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp044_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp046_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp113_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp116_meth%`, na.rm=T)))
#      run=FALSE
#      win_start=win_end
#    }
#    else {
#      All_clones_windowed = rbind(All_clones_windowed, data.frame(chromosome, win_start, win_end,
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp001_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp005_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp044_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp046_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp113_meth%`, na.rm=T),
#      mean(All_cloneschr[(win_start <= All_cloneschr$Position & All_cloneschr$Position <= win_end),]$`SwAsp116_meth%`, na.rm=T)))
#      win_start=win_end
#    }
#  }
#}

#names(All_clones_windowed)<-c("Chromosome","Window_start", "Window_end", "SwAsp001_meth%", "SwAsp005_meth%",
                              "SwAsp044_meth%", "SwAsp046_meth%", "SwAsp113_meth%", "SwAsp116_meth%")

#Calculate mean and standard deviation per window

#All_clones_windowed$Mean<-apply(All_clones_windowed[,4:9],1,mean, na.rm=T)
#All_clones_windowed$STD_dev<-apply(All_clones_windowed[,4:9],1,sd, na.rm=T)

#write.table(All_clones_windowed, "/path/to/filtered/covs/All_chromosomes.windowed.bismark.cov", sep = "\t", quote = FALSE, row.names = FALSE)

#Normalize windows

#All_clones_windowed$`SwAsp001_meth%` = All_clones_windowed$`SwAsp001_meth%`/mean(All_cloneschr$`SwAsp001_meth%`, na.rm=T)
#All_clones_windowed$`SwAsp005_meth%` = All_clones_windowed$`SwAsp005_meth%`/mean(All_cloneschr$`SwAsp005_meth%`, na.rm=T)
#All_clones_windowed$`SwAsp044_meth%` = All_clones_windowed$`SwAsp044_meth%`/mean(All_cloneschr$`SwAsp044_meth%`, na.rm=T)
#All_clones_windowed$`SwAsp046_meth%` = All_clones_windowed$`SwAsp046_meth%`/mean(All_cloneschr$`SwAsp046_meth%`, na.rm=T)
#All_clones_windowed$`SwAsp113_meth%` = All_clones_windowed$`SwAsp113_meth%`/mean(All_cloneschr$`SwAsp113_meth%`, na.rm=T)
#All_clones_windowed$`SwAsp116_meth%` = All_clones_windowed$`SwAsp116_meth%`/mean(All_cloneschr$`SwAsp116_meth%`, na.rm=T)

#save.image("/path/to/filtered/covs/All_clones_window.RData")

#Calculate mean and standard deviation per normalized window

#All_clones_windowed_normalized$Mean<-apply(All_clones_windowed[,4:9],1,mean, na.rm=T)
#All_clones_windowed_normalized$STD_dev<-apply(All_clones_windowed[,4:9],1,sd, na.rm=T)

#Plot to see how windowed methylation works
#for (chromosome in unique(All_clones_windowed$Chromosome)) {
#  All_clones_windowed_chr = All_clones_windowed[All_clones_windowed$Chromosome==chromosome,]
#  plot(All_clones_windowed_chr$Window_start, All_clones_windowed_chr$Mean, main = chromosome, type="l")
#  abline(1,0)

#  plot(All_clones_windowed_chr$Window_start, All_clones_windowed_chr$STD_dev, main = chromosome, type="l")
#}

#Write out table with windowed and normalized data
#write.table(All_clones_windowed, "/path/to/filtered/covs/All_chromosomes.windowed.normalized.bismark.cov", sep = "\t", quote = FALSE, row.names = FALSE)
