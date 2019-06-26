read.table("probes_to_trichocarpa.txt",head=F)->probes_tricho
read.table("consensus_map.txt",head=T)->consensus

names(probes_tricho)<-c("probe","match","pident","length","mismatch","gapopen","pstart","pend","sstart","send","evalue","bitscore")

consensus$scaffold<-as.factor(t(as.data.frame(strsplit(as.character(consensus$marker),":")))[,1])
consensus$position<-as.integer(t(as.data.frame(strsplit(as.character(consensus$marker),":")))[,2])

probes_tricho$tremula_scaffold<-as.factor(t(as.data.frame(strsplit(as.character(probes_tricho$probe),":")))[,1])
probes_tricho$probe_pos<-as.factor(t(as.data.frame(strsplit(as.character(probes_tricho$probe),":")))[,2])
probes_tricho$probe_start<-as.integer(t(as.data.frame(strsplit(as.character(probes_tricho$probe_pos),"-")))[,1])
probes_tricho$probe_end<-as.integer(t(as.data.frame(strsplit(as.character(probes_tricho$probe_pos),"-")))[,2])


###AGAINST TRICHOCARPA v3.0
trichocarpa_map<-NULL
for(i in 1:dim(consensus)[1]){
  df<-probes_tricho[probes_tricho$tremula_scaffold == as.character(consensus$scaffold[i]) & probes_tricho$probe_start <= consensus$position[i] & probes_tricho$probe_end >= consensus$position[i],]
  trichocarpa_map<-rbind(trichocarpa_map,merge(df[df$bitscore==max(df$bitscore),],consensus[i,],by.x="tremula_scaffold", by.y="scaffold"))
} 

LG_Chr_correspondance<-table(trichocarpa_map$LG,trichocarpa_map$match)[,1:19]


###Create synteny map between tremula genetic map and trichocarpa physical map

pdf("tremula genetic distance vs. trichocarpa physical distance.pdf")
par(mfrow=c(19,19))
par(mar=c(0,0,0,0))
par(omi=c(1,1,0,0))
par(xpd=TRUE)
mismatches<-NULL
for(i in 1:19){
  for(j in 1:19){
    if(i %in% 1:9){
      tmp<-subset(trichocarpa_map,LG==j & match==paste("Chr0",i,sep=""))
    }else{
      tmp<-subset(trichocarpa_map,LG==j & match==paste("Chr",i,sep=""))
    }
    tmp$trich_pos<-tmp$sstart + (tmp$position - tmp$probe_start)
    #if(dim(tmp)[1]>150){
    #  print(cor(tmp$consensus,tmp$trich_pos,meth="kend"))
    #}else{
    #  mismatches<-c(mismatches,dim(tmp)[1])
   # }
    if(dim(tmp)[1]>0){
      col<-ifelse(dim(tmp)[1]>150,"black","grey65")
      plot(trich_pos~consensus,tmp,xlab="",ylab="",col=col,pch=19,xaxt="n",yaxt="n",cex=0.7)
    }else{
      plot(c(0,0),c(0,0),type="n",xlab="",ylab="",col=col,pch=19,xaxt="n",yaxt="n",cex=0.7)
    }
  }
}
mtext(c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11","LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19"),side=1,line=1, outer=TRUE, at=seq(0.025,1,0.053),las=2)
mtext(c("Chr19","Chr18","Chr17","Chr16","Chr15","Chr14","Chr13","Chr12","Chr11","Chr10","Chr09","Chr08","Chr07","Chr06","Chr05","Chr04","Chr03","Chr02","Chr01"),side=2,line=1, outer=TRUE, at=seq(0.03,1,0.052),las=1)
mtext(expression(paste(italic("P. tremula"),", Consensus map")),side=1, line=5,outer=T,at=0.5)
mtext(expression(paste(italic("P. trichocarpa"),", assembly v3.0, hardmasked")),side=2, line=6,outer=T,at=0.5)
dev.off()
rm(tmp)