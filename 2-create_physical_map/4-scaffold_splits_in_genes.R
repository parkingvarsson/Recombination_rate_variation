read.table("Potra-v1.1-annotation.nocom.onecopy.deduplicated.sorted.with_weird_gaps.gff",head=F)->strange_gap_mrna
read.table("scaffold_length.txt",head=F)->Scaffold_length
read.table("breakpoints_final.bed",head=F)->breakpoints

Scaffold_length<-Scaffold_length[Scaffold_length$V1 %in% unique(as.character(strange_gap_mrna$V1)),]
breakpoints<-breakpoints[breakpoints$V1%in% unique(as.character(strange_gap_mrna$V1)),]

library(gdata)
Scaffold_length<-drop.levels(Scaffold_length)
breakpoints<-drop.levels(breakpoints)

pdf("strangly_split_scaffolds.pdf")
par(mfrow=c(3,2))
for(i in 1:length(unique(strange_gap_mrna$V1))){
  scaf<-unique(as.character(strange_gap_mrna$V1))[i]
  breaks<-breakpoints[breakpoints$V1 == scaf,]
  mrna_plus<-strange_gap_mrna[strange_gap_mrna$V1==scaf & strange_gap_mrna$V7=="+",]
  mrna_minus<-strange_gap_mrna[strange_gap_mrna$V1==scaf& strange_gap_mrna$V7=="-",]
  length<-Scaffold_length$V2[Scaffold_length$V1==scaf]
  plot(NULL,ylim=c(0,2),xlim=c(0,length),main=scaf,ylab="",xlab="Genomic position, bp",yaxt="n")
  for(l in 1:dim(breaks)[1]){
    rect(breaks$V2[l],0,breaks$V3[l],2,col=rgb(red = 0,green = 0, blue = 1, alpha=0.5),border=rgb(red = 0,green = 0, blue = 1, alpha=0.5))
  }
  if(dim(mrna_plus)[1]>0){
    for(j in 1:length(mrna_plus$V4[mrna_plus$V3=="mRNA"])){
      lines(x=c(mrna_plus$V4[mrna_plus$V3=="mRNA"][j],mrna_plus$V5[mrna_plus$V3=="mRNA"][j]),y=c(1.5,1.5))
    }
    for(m in 1:length(mrna_plus$V4[mrna_plus$V3=="exon"])){
      rect(mrna_plus$V4[mrna_plus$V3=="exon"][m],1.3,mrna_plus$V5[mrna_plus$V3=="exon"][m],1.7,col=rgb(red = 1,green = 0, blue = 0, alpha=0.5))
    }
  }
  if(dim(mrna_minus)[1]>0){
    for(k in 1:length(mrna_minus$V4[mrna_minus$V3=="mRNA"])){
      lines(x=c(mrna_minus$V4[mrna_minus$V3=="mRNA"][k],mrna_minus$V5[mrna_minus$V3=="mRNA"][k]),y=c(0.5,0.5))
    }
    for(n in 1:length(mrna_minus$V4[mrna_minus$V3=="exon"])){
      rect(mrna_minus$V4[mrna_minus$V3=="exon"][n],0.3,mrna_minus$V5[mrna_minus$V3=="exon"][n],0.7,col=rgb(red = 0,green = 1, blue = 0, alpha=0.5))
    }
  }  
 }  
dev.off()
rm(scaf,breaks,mrna_minus,mrna_plus,length)