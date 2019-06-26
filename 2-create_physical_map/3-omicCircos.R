######################################
##Make a circular plot of linkage groups and connect markers from the same scaffolds
######################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("OmicCircos")

read.table("Consensus_with_bins.txt",head=T)->maps
maps$scaffold<-NULL
for(i in 1:dim(maps)[1]){
  maps$scaffold[i]<-strsplit(as.character(maps$marker[i]),split=":")[[1]][1]
} 
maps$scaffold<-as.factor(maps$scaffold)

maps<-maps[order(maps$Chr,maps$consensus),]

library (OmicCircos)
options (stringsAsFactors=FALSE) ;
set.seed(1234) ;

##Segment frame
cons.name <- paste("Chr",maps$Chr,sep=" ")
cons.start <- as.character(round(maps$consensus),digit=3)
cons.end <- as.character(round(maps$consensus),digit=3)
cons.f <- data.frame(seg.name=cons.name, seg.start=cons.start,seg.end=cons.end,the.v=runif(length(cons.end)),Note=as.character(maps$marker))
db <- segAnglePo(cons.f,seg=unique(cons.name),angle.start = 3,angle.end = 357);

##Segment mapping
cons.v<-maps[,c(4,2)]
names(cons.v)<-c("seg.name","seg.pos")
cons.v$seg.pos<-as.character(cons.v$seg.pos)
cons.v$seg.name<-paste("Chr",cons.v$seg.name,sep=" ")

###Segment links
cons.links <- NULL
cons.links$chr1 <-NULL
cons.links$pos1 <-NULL
cons.links$marker1 <-NULL
cons.links$chr2 <-NULL
cons.links$pos2 <-NULL
cons.links$marker2 <-NULL
for (i in 1:length(levels(maps$scaffold))){
  tmp <- maps[maps$scaffold == levels(maps$scaffold)[i],]
  if (dim(tmp)[1] > 1){
    for (j in 2:dim(tmp)[1]){
      cons.links$chr1 <- c(cons.links$chr1,as.character(paste("Chr",tmp$Chr[j-1],sep=" ")))
      cons.links$chr2 <- c(cons.links$chr2,as.character(paste( "Chr", tmp$Chr[j],sep=" ")))
      cons.links$pos1 <- c(cons.links$pos1,as.character(tmp$consensus[j-1]))
      cons.links$pos2 <- c(cons.links$pos2,as.character(tmp$consensus[j]))
      cons.links$marker1 <- c(cons.links$marker1,as.character(tmp$marker[j-1]))
      cons.links$marker2 <- c(cons.links$marker2,as.character(tmp$marker[j]))
    }
  }
}

cons.links <- as.data.frame(cons.links)
cons.links <- cons.links[,c(1,3,5,2,4,6)]

intra_cons.links <- cons.links[cons.links$chr1 == cons.links$chr2,]
inter_cons.links <- cons.links[cons.links$chr1 != cons.links$chr2,]

col.intra.links<-NULL
for ( i in 1:dim(intra_cons.links)[1]){
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) >=20)){
    col.intra.links[i] <- "red"
  }
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) < 20)){
    col.intra.links[i] <- "grey90"
  }
}


col.inter.links<-NULL
tmp<-t(as.data.frame(strsplit(inter_cons.links$marker1,":"))[1,])
tmp2<-t(as.data.frame(strsplit(inter_cons.links$marker2,":"))[1,])
tmp<-table(tmp)
tmp2<-table(tmp2)
tmp<-unique(names(c(tmp[tmp>1],tmp2[tmp2>1])))
for (i in 1:dim(inter_cons.links)[1]){
  if (strsplit(inter_cons.links$marker1[i],":")[[1]][1] %in% tmp){
    col.inter.links[i]<-"dark blue"
  } 
  else{
    col.inter.links[i]<-"orange"  
  } 
}

inter_cons.links2<-inter_cons.links[which(col.inter.links == "dark blue"),]
intra_cons.links2<-intra_cons.links[which(col.intra.links== "red"),]



pdf("Circosplot_with_split_Scaffolds.pdf")
par(mar=c(2,2,2,2),mfrow=c(1,1))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
circos(R=350, cir=db, type="chr",col="grey50", print.chr.lab=TRUE, W=4, scale=FALSE)
circos(R=310, cir=db, mapping=cons.v, type="b3",col="black", W=40 ,lwd=1,B=TRUE,scale=FALSE)
circos(R=305, cir=db, mapping=intra_cons.links, type="link2",B=TRUE, col=col.intra.links)
circos(R=305, cir=db, mapping=intra_cons.links2, type="link2", col="red")
circos(R=185, cir=db, mapping=inter_cons.links, type="link",B=TRUE, col=col.inter.links)
circos(R=185, cir=db, mapping=inter_cons.links2, type="link", col="dark blue")
text(400,730,"A")
text(400,650,"B")
text(400,550,"C")
dev.off()
