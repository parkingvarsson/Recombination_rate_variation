read.table("d1.10-Maps.txt",head=F)->d1.10
read.table("d2.15-Maps.txt",head=F)->d2.15
read.table("MultiMarkerBins_with_MapPositions.txt",head=T)->bin_pos

tmp<-bin_pos[!is.na(bin_pos$D1.10_LG) & !(bin_pos$markers %in% d1.10$V2),c(4,1,5)]
names(tmp)<-c("V1","V2","V3")
d1.10_with_bins<-rbind(d1.10,tmp)
tmp<-bin_pos[!is.na(bin_pos$D2.15_LG)& !(bin_pos$markers %in% d2.15$V2),c(6,1,7)]
names(tmp)<-c("V1","V2","V3")
d2.15_with_bins<-rbind(d2.15,tmp)
rm(tmp) 

library(gdata)
d1.10_with_bins<-drop.levels(d1.10_with_bins)
d2.15_with_bins<-drop.levels(d2.15_with_bins)

d1.10_with_bins<-d1.10_with_bins[order(d1.10_with_bins$V1,d1.10_with_bins$V3),]
d2.15_with_bins<-d2.15_with_bins[order(d2.15_with_bins$V1,d2.15_with_bins$V3),]

d1.10_with_bins$V1<-as.factor(d1.10_with_bins$V1)
d2.15_with_bins$V1<-as.factor(d2.15_with_bins$V1)

d1.10_with_bins$V4<-as.factor(t(as.data.frame(strsplit(as.character(d1.10_with_bins$V2),":")))[,1])
d1.10_with_bins$V5<-as.factor(t(as.data.frame(strsplit(as.character(d1.10_with_bins$V2),":")))[,2])

d2.15_with_bins$V4<-as.factor(t(as.data.frame(strsplit(as.character(d2.15_with_bins$V2),":")))[,1])
d2.15_with_bins$V5<-as.factor(t(as.data.frame(strsplit(as.character(d2.15_with_bins$V2),":")))[,2])

names(d1.10_with_bins)<-c("LG","Marker","genetic position","Scaffold ID", "scaffold position")
names(d2.15_with_bins)<-c("LG","Marker","genetic position","Scaffold ID", "scaffold position")


length(unique(c(levels(d1.10_with_bins$`Scaffold ID`),levels(d2.15_with_bins$`Scaffold ID`)))) ### 4351 scaffolds combined
length(unique(c(levels(d1.10_with_bins$Marker),levels(d2.15_with_bins$Marker)))) ### 19519 markers combined

LGs_2_chr<-as.data.frame(cbind(c(seq(1,19,1)),c(4,9,2,10,14,1,5,3,19,15,13,8,18,6,17,11,12,7,16)))
names(LGs_2_chr)<-c("LG","Chr")

d1.10_with_bins$Chr<-as.factor(unlist(lapply(as.numeric(d1.10_with_bins$LG),function(x) LGs_2_chr$Chr[LGs_2_chr$LG==x])))
d2.15_with_bins$Chr<-as.factor(unlist(lapply(as.numeric(d2.15_with_bins$LG),function(x) LGs_2_chr$Chr[LGs_2_chr$LG==x])))

d1.10_with_bins<-d1.10_with_bins[order(d1.10_with_bins$Chr,d1.10_with_bins$`genetic position`),c(6,1,2,3,4,5)]
d2.15_with_bins<-d2.15_with_bins[order(d2.15_with_bins$Chr,d2.15_with_bins$`genetic position`),c(6,1,2,3,4,5)]

maps_combined<-merge(d1.10_with_bins,d2.15_with_bins,by=c("Marker","Chr","LG","Scaffold ID","scaffold position"),all = T)
maps_combined<-maps_combined[order(maps_combined$Chr,maps_combined$`genetic position.x`),]
maps_combined$`scaffold position`<-as.integer(as.character(maps_combined$`scaffold position`))

split_scafs<-NULL
for(i in 1:18){
  for(j in (i+1):19){
    split_scafs<-c(split_scafs,as.character(unique(maps_combined$`Scaffold ID`[maps_combined$Chr == i][maps_combined$`Scaffold ID`[maps_combined$Chr == i] %in% maps_combined$`Scaffold ID`[maps_combined$Chr == j]])))
  }
}
split_scafs<-unique(split_scafs)

Inter_split_map<-maps_combined[maps_combined$`Scaffold ID` %in% split_scafs,]
Inter_split_map<-drop.levels(Inter_split_map)

Inter_scaf_gap_bed<-NULL
for(i in 1:length(levels(Inter_split_map$`Scaffold ID`))){
  tmp<-Inter_split_map[Inter_split_map$`Scaffold ID` == levels(Inter_split_map$`Scaffold ID`)[i],]
  tmp<-tmp[order(tmp$`scaffold position`),]
  for(j in 2:dim(tmp)[1]){
    if(tmp$LG[j-1] != tmp$LG[j]){
      Inter_scaf_gap_bed<-rbind(Inter_scaf_gap_bed,c(as.character(tmp$`Scaffold ID`[j]),tmp$`scaffold position`[j-1],tmp$`scaffold position`[j]))
    }
  }
}
Inter_scaf_gap_bed<-as.data.frame(Inter_scaf_gap_bed)
Inter_scaf_gap_bed$V2<-as.integer(as.character(Inter_scaf_gap_bed$V2))
Inter_scaf_gap_bed$V3<-as.integer(as.character(Inter_scaf_gap_bed$V3))


Nr_of_LG<-NULL
for (i in 1:length(split_scafs)){
  Nr_of_LG<-c(Nr_of_LG, length(unique(maps_combined$Chr[maps_combined$`Scaffold ID`==split_scafs[i]])))
}



read.table("Scaffold_length.txt",head=F)->genome_scaf_length
read.table("Potra01-geneOnly.gff3",head=F)->genes
genes[genes$V1 %in% split_scafs,]->genes
read.table("v1.1_Potra01-genome_gaps.bed",head=F)->assembly_gaps
assembly_gaps[assembly_gaps$V1 %in% split_scafs,]->assembly_gaps
assembly_gaps

split_pos_Inter_split_scaf<-NULL
no_Inter_split_pos_availible<-NULL
for(i in 1:length(split_scafs)){
  map_split_pos<-Inter_scaf_gap_bed[Inter_scaf_gap_bed$V1==split_scafs[i],]
  for(j in 1:dim(map_split_pos)[1]){
    Intergenic_assembly_splits<-as.data.frame(NULL)
    tmp<-NULL
    tmp1<-NULL
    assembly_splits<-assembly_gaps[assembly_gaps$V1==split_scafs[i] & assembly_gaps$V2 > map_split_pos$V2[j] & assembly_gaps$V3 < map_split_pos$V3[j],]
    assembly_splits$V4<-assembly_splits$V3-assembly_splits$V2
    gene_pos<-genes[genes$V1 == split_scafs[i] & genes$V5 > map_split_pos$V2[j] & genes$V4 < map_split_pos$V3[j],]
    
    if(dim(assembly_splits)[1] != 0 & dim(gene_pos)[1]!=0){
      tmp<-matrix(NA,nrow = dim(assembly_splits)[1],ncol = dim(gene_pos)[1],byrow = F)
      for(y in 1:dim(assembly_splits)[1]){
        for(z in 1:dim(gene_pos)[1]){
          tmp[y,z]<-assembly_splits$V2[y] < gene_pos$V4[z] |assembly_splits$V2[y] > gene_pos$V5[z] & assembly_splits$V3[y] < gene_pos$V4[z] |assembly_splits$V3[y] > gene_pos$V5[z]
        }
      }
      tmp1<-apply(as.data.frame(tmp),1,function(y) sum(y=="TRUE"))
      Intergenic_assembly_splits<-assembly_splits[which(tmp1==dim(as.data.frame(tmp))[2]),]
      Intergenic_assembly_splits<-Intergenic_assembly_splits[!is.na(Intergenic_assembly_splits$V2),]
    }
    if(dim(assembly_splits)[1] == 0){
      no_Inter_split_pos_availible<-c(no_Inter_split_pos_availible,split_scafs[i])
      ### if no assembly gaps are availible and no genes are in the split region, split the scaffold in the middle of the split region
      if(dim(gene_pos)[1] == 0){
        split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(as.numeric(map_split_pos[j,2:3]))),ceiling(median(as.numeric(map_split_pos[j,2:3])))+1,1))))
      }else{
        ### If genes are in the split region, find a region outside the gene region and make the split in the middle of this region
        if(dim(gene_pos)[1]>1){
          if(gene_pos$V4[2]>gene_pos$V5[1]){
            split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(c(as.numeric(gene_pos$V5[1]),as.numeric(gene_pos$V4[2])))),ceiling(median(c(as.numeric(gene_pos$V5[1]),as.numeric(gene_pos$V4[2]))))+1,1))))
          }else{
            if(min(gene_pos$V4) > map_split_pos$V2[j]){
            split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(c(as.numeric(map_split_pos$V2[j]),as.numeric(min(gene_pos$V4))))),ceiling(median(c(as.numeric(map_split_pos$V2[j]),as.numeric(min(gene_pos$V4)))))+1,1))))
            }else{
              if(max(gene_pos$V5) < map_split_pos$V3[j]){
                split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(c(as.numeric(map_split_pos$V3[j]),as.numeric(max(gene_pos$V5))))),ceiling(median(c(as.numeric(map_split_pos$V3[j]),as.numeric(max(gene_pos$V5)))))+1,1))))
              }else{
                ###If a gene covers the entire split region, split the scaffold in the middle of the split region
                split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(as.numeric(map_split_pos[j,2:3]))),ceiling(median(as.numeric(map_split_pos[j,2:3])))+1,1))))
              }
            }
          }
        }else{
          if(gene_pos$V4 > map_split_pos$V2[j]){
            split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(c(as.numeric(map_split_pos$V2[j]),as.numeric(gene_pos$V4[2])))),ceiling(median(c(as.numeric(map_split_pos$V2[j]),as.numeric(gene_pos$V4[2]))))+1,1))))
          }else{
            if(gene_pos$V5 < map_split_pos$V3[j]){
              split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(c(as.numeric(map_split_pos$V3[j]),as.numeric(gene_pos$V5[1])))),ceiling(median(c(as.numeric(map_split_pos$V3[j]),as.numeric(gene_pos$V5[1]))))+1,1))))
            }else{
              ###If a gene covers the entire split region, split the scaffold in the middle of the split region
              split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,as.data.frame(t(c(as.character(map_split_pos$V1[j]),ceiling(median(as.numeric(map_split_pos[j,2:3]))),ceiling(median(as.numeric(map_split_pos[j,2:3])))+1,1))))
            }
          }  
        }
      }  
    }else{
      if(dim(Intergenic_assembly_splits)[1]!=0){
        split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,Intergenic_assembly_splits[Intergenic_assembly_splits$V4 == max(Intergenic_assembly_splits$V4),])
      }else{
        split_pos_Inter_split_scaf<-rbind(split_pos_Inter_split_scaf,assembly_splits[assembly_splits$V4 == max(assembly_splits$V4),])
      } 
    }
  }
}

no_Inter_split_pos_availible<-unique(no_Inter_split_pos_availible)

input<-sort(split_scafs)
#input<-as.character(split_pos_Inter_split_scaf$V1)
#input<-no_Inter_split_pos_availible
#input<-no_Inter_split_pos_availible[no_Inter_split_pos_availible %in% split_pos_Inter_split_scaf$V1]
#input<-no_Inter_split_pos_availible[!(no_Inter_split_pos_availible %in% split_pos_Inter_split_scaf$V1)]



pdf("inter_split_scaffolds_split_decisions.pdf")
par(mfrow=c(3,3))
plot.new()
legend("center",pch=c(19,NA,NA,NA,NA),lwd=c(NA,1,1,1,1),col=c("purple","black","green","blue","red"),c("Probe-marker","Gene","Assembly gaps","Gap split", "Artificial split"),bty="n")
for(i in 1:length(input)){
  plot(Inter_split_map$`scaffold position`[Inter_split_map$`Scaffold ID` == input[i]],as.integer(Inter_split_map$Chr[Inter_split_map$`Scaffold ID` == input[i]]),las=1,pch=19,col="purple",cex=0.5,xlab="Genomic position (bp)",ylab="Chromosome",ylim=c(0,19),main=as.character(input[i]),xlim=c(0,genome_scaf_length$V2[genome_scaf_length$V1==input[i]]))
  #points(as.integer(as.character(d2.15_with_bins$`scaffold position`[d2.15_with_bins$`Scaffold ID` == input[i]])),as.integer(d2.15_with_bins$Chr[d2.15_with_bins$`Scaffold ID` == input[i]])+0.2,pch=19,col="blue",cex=0.5)
  df<-genes[genes$V1 == input[i],]
  for(j in 1:length(df$V4)){
    lines(c(as.numeric(df$V4[j]),as.numeric(df$V5[j])),c(0,0),col="black")
  }
  tmp<-assembly_gaps[assembly_gaps$V1==input[i],]
  tmp$V4<-tmp$V3-tmp$V2
  tmp1<-split_pos_Inter_split_scaf[split_pos_Inter_split_scaf$V1==input[i],]
  for(k in 1:dim(tmp)[1]){
    if (dim(tmp)[1] != 0 ){
      rect(tmp$V2[k],0.5,tmp$V3[k],19,density=NA, col=rgb(red=0, green=0.7, blue=0,alpha=0.5))
    }
  }
  for(l in 1:dim(tmp1)[1]){
    if(tmp1$V4[l]==1){
      rect(as.numeric(tmp1$V2[l]),0.5,as.numeric(tmp1$V3[l]),19,density=NA, col=rgb(red=1, green=0, blue=0,alpha=0.5))
    }else{
      rect(as.numeric(tmp1$V2[l]),0.5,as.numeric(tmp1$V3[l]),19,density=NA, col=rgb(red=0, green=0, blue=1,alpha=0.5))
    }
  }
  
  abline(h=seq(1,19,1),lty=3,lwd=0.3,col="grey50")
}
dev.off()


### Intra chromosomes split scaffolds
par(mfrow=c(3,3))
Intra_split_scaffolds<-NULL
for(i in 1:length(levels(maps_combined$`Scaffold ID`))){
  scaffold_map<-maps_combined[maps_combined$`Scaffold ID`==levels(maps_combined$`Scaffold ID`)[i],]
  for(j in 1:length(unique(scaffold_map$Chr))){
    scaffold_Chr_map<-scaffold_map[scaffold_map$Chr==unique(scaffold_map$Chr)[j],]
    if(max(scaffold_Chr_map$`genetic position.x`,na.rm=T) - min(scaffold_Chr_map$`genetic position.x`,na.rm=T) > 20 |max(scaffold_Chr_map$`genetic position.y`,na.rm=T) - min(scaffold_Chr_map$`genetic position.y`,na.rm=T) > 20){
      Intra_split_scaffolds<-c(Intra_split_scaffolds,as.character(unique(scaffold_Chr_map$`Scaffold ID`)))
      #plot(scaffold_Chr_map$`scaffold position`,scaffold_Chr_map$`genetic position.x`,ylim=c(0,max(as.numeric(maps_combined$`genetic position.x`[maps_combined$Chr==unique(as.character(scaffold_Chr_map$Chr))]),na.rm=T)),xlim=c(0,genome_scaf_length$V2[genome_scaf_length$V1==unique(as.character(scaffold_Chr_map$`Scaffold ID`))]),las=1,xlab="Genomic position (bp)",ylab="Genetic distance (cM)",main=paste(unique(as.character(scaffold_Chr_map$`Scaffold ID`)),", Chr",unique(scaffold_Chr_map$Chr,sep="")),pch=19,cex=0.7,col="red")
      #points(scaffold_Chr_map$`scaffold position`,scaffold_Chr_map$`genetic position.y`,col="blue",cex=0.7,pch=19)
    }
  }
}

Intra_split_scaffold_map<-NULL
for(i in 1:length(Intra_split_scaffolds)){
  Intra_split_scaffold_map<-rbind(Intra_split_scaffold_map,maps_combined[maps_combined$`Scaffold ID`==Intra_split_scaffolds[i],][order(maps_combined[maps_combined$`Scaffold ID`==Intra_split_scaffolds[i],]$`scaffold position`),])
}

Intra_split_scaffold_map[Intra_split_scaffold_map$`Scaffold ID` %in% split_scafs,]
Intra_split_scaffold_map[!(Intra_split_scaffold_map$`Scaffold ID` %in% split_scafs),]


###Intra_split_regions.txt got manually evaluated. Therefore read in the table
read.table("Intra_split_regions.txt",head=T)-> Intra_split_regions

no_Intra_split_pos_availible<-NULL
split_pos_Intra_split_scaf<-NULL
for (i in 1:dim(Intra_split_regions)[1]){
  map_split_pos<-Intra_split_regions[i,]
  assembly_splits<-assembly_gaps[assembly_gaps$V1==as.character(map_split_pos$Scaffold) & assembly_gaps$V2 > map_split_pos$Split_start & assembly_gaps$V3 < map_split_pos$Split_end,]
  assembly_splits$V4<-assembly_splits$V3 - assembly_splits$V2
  gene_pos<-genes[genes$V1==as.character(map_split_pos$Scaffold) & genes$V5 > map_split_pos$Split_start & genes$V4 < map_split_pos$Split_end,]
  tmp<-NULL
  tmp1<-NULL
  Intergenic_assembly_splits<-as.data.frame(NULL)
  if(dim(assembly_splits)[1] > 0 & dim(gene_pos)[1] > 0){
    tmp<-matrix(NA,nrow = dim(assembly_splits)[1],ncol = dim(gene_pos)[1],byrow = F)
    for(y in 1:dim(assembly_splits)[1]){
      for(z in 1:dim(gene_pos)[1]){
        tmp[y,z]<-assembly_splits$V2[y] < gene_pos$V4[z] |assembly_splits$V2[y] > gene_pos$V5[z] & assembly_splits$V3[y] < gene_pos$V4[z] |assembly_splits$V3[y] > gene_pos$V5[z]
      }
    }
    tmp1<-apply(as.data.frame(tmp),1,function(y) sum(y=="TRUE"))
    Intergenic_assembly_splits<-assembly_splits[which(tmp1==dim(as.data.frame(tmp))[2]),]
    Intergenic_assembly_splits<-Intergenic_assembly_splits[!is.na(Intergenic_assembly_splits$V2),]
  }
  if(dim(assembly_splits)[1] > 0 & dim(gene_pos)[1] == 0){
    split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,assembly_splits[assembly_splits$V4 == max(assembly_splits$V4),])
  }
  if(dim(assembly_splits)[1] == 0){
    no_Intra_split_pos_availible<-c(no_Intra_split_pos_availible,as.character(map_split_pos$Scaffold))
    if(dim(gene_pos)[1]==0){
      split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(as.numeric(map_split_pos[2:3]))),ceiling(median(as.numeric(map_split_pos[2:3])))+1,1))))
    }else{
      ### If genes are in the split region, find a region outside the gene region and make the split in the middle of this region
      if(dim(gene_pos)[1] > 1){
        if(gene_pos$V4[2] > gene_pos$V5[1]){
          split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(c(as.numeric(gene_pos$V5[1]),as.numeric(gene_pos$V4[2])))),ceiling(median(c(as.numeric(gene_pos$V5[1]),as.numeric(gene_pos$V4[2]))))+1,1))))
        }else{
          if(min(gene_pos$V4) > map_split_pos$Split_start){
            split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(c(as.numeric(map_split_pos$Split_start),as.numeric(min(gene_pos$V4))))),ceiling(median(c(as.numeric(map_split_pos$Split_start),as.numeric(min(gene_pos$V4)))))+1,1))))
          }else{
            if(max(gene_pos$V5) < map_split_pos$Split_end){
              split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(c(as.numeric(map_split_pos$Split_end),as.numeric(max(gene_pos$V5))))),ceiling(median(c(as.numeric(map_split_pos$Split_end),as.numeric(max(gene_pos$V5)))))+1,1))))
            }else{
              ###If a gene covers the entire split region, split the scaffold in the middle of the split region
              split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(as.numeric(map_split_pos[2:3]))),ceiling(median(as.numeric(map_split_pos[2:3])))+1,1))))
            }
          }
        }
      }else{
        if(gene_pos$V4 > map_split_pos$Split_start){
          split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(c(as.numeric(map_split_pos$Split_start),as.numeric(gene_pos$V4[2])))),ceiling(median(c(as.numeric(map_split_pos$Split_start),as.numeric(gene_pos$V4[2]))))+1,1))))
        }else{
          if(gene_pos$V5 < map_split_pos$Split_end){
            split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(c(as.numeric(map_split_pos$Split_end),as.numeric(gene_pos$V5[1])))),ceiling(median(c(as.numeric(map_split_pos$Split_end),as.numeric(gene_pos$V5[1]))))+1,1))))
          }else{
            ###If a gene covers the entire split region, split the scaffold in the middle of the split region
            split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,as.data.frame(t(c(as.character(map_split_pos$Scaffold),ceiling(median(as.numeric(map_split_pos[2:3]))),ceiling(median(as.numeric(map_split_pos[2:3])))+1,1))))
          }
        }  
      }
    }
  }else{
    if(dim(Intergenic_assembly_splits)[1] > 0){
      split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,Intergenic_assembly_splits[Intergenic_assembly_splits$V4 == max(Intergenic_assembly_splits$V4),])
    }else{
      split_pos_Intra_split_scaf<-rbind(split_pos_Intra_split_scaf,assembly_splits[assembly_splits$V4 == max(assembly_splits$V4),])
    }
  }  
}

input<-Intra_split_scaffolds[Intra_split_scaffolds %in% as.character(split_pos_Intra_split_scaf$V1)]
pdf("intra_split_scaffolds_split_decisions.pdf")
par(mfrow=c(3,3))
plot.new()
legend("center",pch=c(19,19,19,NA,NA,NA,NA),lwd=c(NA,NA,NA,1,1,1,1),col=c("red","blue","grey","black","green","blue","red"),c("Female probe-marker","Male probe-marker","Probe-marker on other Chr","Gene","Assembly gaps","Split gap", "Artificial gap"),bty="n",cex=0.9)
for(i in 1:length(input)){
  tm<-Intra_split_scaffold_map[Intra_split_scaffold_map$`Scaffold ID` == input[i],]
  size<-max(c(tm$`genetic position.x`,tm$`genetic position.y`),na.rm=T)
  if (length(unique(tm$Chr))==1){
    plot(tm$`scaffold position`,tm$`genetic position.x`,las=1,pch=19,col="red",cex=0.5,xlab="Genomic position (bp)",ylab="Genetic distance",main=paste(as.character(input[i])," - Chr",unique(as.character(tm$Chr)),sep=""),ylim=c(-10,size),xlim=c(0,genome_scaf_length$V2[genome_scaf_length$V1==input[i]]))
    points(tm$`scaffold position`,tm$`genetic position.y`,pch=19,col="blue",cex=0.5)
  }else{
    for(m in 1:length(unique(as.character(tm$Chr)))){
      if (max(c(tm$`genetic position.x`[tm$Chr==unique(as.character(tm$Chr))[m]],tm$`genetic position.y`[tm$Chr==unique(as.character(tm$Chr))[m]]),na.rm=T)-min(c(tm$`genetic position.x`[tm$Chr==unique(as.character(tm$Chr))[m]],tm$`genetic position.y`[tm$Chr==unique(as.character(tm$Chr))[m]]),na.rm=T) > 20 ){
        z<-unique(as.character(tm$Chr))[m]
      }
    }
    plot(tm$`scaffold position`[tm$Chr == z],tm$`genetic position.x`[tm$Chr == z],las=1,pch=19,col="red",cex=0.5,xlab="Genomic position (bp)",ylab="Genetic distance",main=paste(as.character(input[i])," - Chr",z,sep=""),ylim=c(-10,size),xlim=c(0,genome_scaf_length$V2[genome_scaf_length$V1==input[i]]))
    points(tm$`scaffold position`[tm$Chr == z],tm$`genetic position.y`[tm$Chr == z],pch=19,col="blue",cex=0.5)
    points(tm$`scaffold position`[tm$Chr != z],tm$`genetic position.x`[tm$Chr != z],pch=19,col="gray50",cex=0.5)
    points(tm$`scaffold position`[tm$Chr != z],tm$`genetic position.y`[tm$Chr != z],pch=19,col="gray50",cex=0.5)
  }
  df<-genes[genes$V1 == input[i],]
  for(j in 1:length(df$V4)){
    lines(c(as.numeric(df$V4[j]),as.numeric(df$V5[j])),c(-10,-10),col="black")
  }
  tmp<-assembly_gaps[assembly_gaps$V1==input[i],]
  tmp$V4<-tmp$V3-tmp$V2
  tmp1<-split_pos_Intra_split_scaf[split_pos_Intra_split_scaf$V1==input[i],]
  for(k in 1:dim(tmp)[1]){
    if (dim(tmp)[1] != 0 ){
      rect(tmp$V2[k],0.5,tmp$V3[k],size,density=NA, col=rgb(red=0, green=0.7, blue=0,alpha=0.5))
    }
  }
  for(l in 1:dim(tmp1)[1]){
    if(tmp1$V4[l]==1){
      rect(as.numeric(tmp1$V2[l]),0,as.numeric(tmp1$V3[l]),size,density=NA, col=rgb(red=1, green=0, blue=0,alpha=0.5))
    }else{
      rect(as.numeric(tmp1$V2[l]),0,as.numeric(tmp1$V3[l]),size,density=NA, col=rgb(red=0, green=0, blue=1,alpha=0.5))
    }
  }
}
dev.off()


split_pos_Intra_split_scaf<-split_pos_Intra_split_scaf[-which(split_pos_Intra_split_scaf$V1=="Potra001073"),]
#####REMOVE following markers from genetic map:
#Potra001073:48123 ##This region has overlapping genes and it is impossible to know what is correct and not
#Potra001073:57711

breakpoints_final<-rbind(split_pos_Inter_split_scaf,split_pos_Intra_split_scaf)
breakpoints_final<-breakpoints_final[order(breakpoints_final$V1,breakpoints_final$V2),]
write.table(breakpoints_final,"breakpoints_final.bed",row.names = F,col.names = F,quote=F,sep="\t")

### ADD BIN MARKERS TO THE CONSENSUS MAP
read.table("Consensus_map.txt",head=T)->Consensus
bin_pos$Consensus<-NA
for(i in 1:length(unique(bin_pos$bin))){
  bin_pos$Consensus[bin_pos$bin==unique(bin_pos$bin)[i]]<-Consensus$consensus[Consensus$marker %in% as.character(bin_pos$markers[bin_pos$bin==unique(bin_pos$bin)[i]])]
}
bin_pos$LG<-NA
for(i in 1:dim(bin_pos)[1]){
  bin_pos$LG[i]<-unique(c(bin_pos$D1.10_LG[i],bin_pos$D2.15_LG[i])[!is.na(c(bin_pos$D1.10_LG[i],bin_pos$D2.15_LG[i]))])
}
Consensus_with_bins<-merge(Consensus,bin_pos[,c(1,8,9)],by.x=c("marker","consensus","LG"),by.y=c("markers","Consensus","LG"),all = T)
Consensus_with_bins<-Consensus_with_bins[order(Consensus_with_bins$LG,Consensus_with_bins$consensus),]
Consensus_with_bins$Chr<-as.factor(unlist(lapply(as.numeric(Consensus_with_bins$LG),function(x) LGs_2_chr$Chr[LGs_2_chr$LG==x])))

write.table(Consensus_with_bins,"Consensus_with_bins.txt",quote=F,col.names = T,row.names = F,sep="\t")