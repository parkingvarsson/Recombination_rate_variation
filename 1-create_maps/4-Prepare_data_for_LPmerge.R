read.table("d1.10-Maps.txt",head=F)->d1.10
read.table("d2.15-Maps.txt",head=F)->d2.15

names(d1.10)<-c("LG","Marker","Size")
names(d2.15)<-c("LG","Marker","Size")


###Control map orientations
plotCompareMap<-function(d1.10,d2.15,LG,max_length){
  
  par(las=1)
  plot(x=c(1,2),y=c(0,max_length),xlim=c(0.90,2.1),type="n",yaxt="n",xaxt="n",xlab="", ylab=expression("Genetic Distance (cM)"),main=paste("LG",LG))
  axis(side=2,at=seq(0,max_length,50),labels=seq(max_length,0,-50))
  axis(side=1,at=c(1,2),labels=c("d1.10","d2.15"),cex=0.7)
  
  lines(x=c(1,1),y=c(max_length-max(d1.10$Size[d1.10$LG==LG],na.rm=T),max_length))
  lines(x=c(2,2),y=c(max_length-max(d2.15$Size[d2.15$LG==LG],na.rm=T),max_length))
  
  
  for(i in 1:dim(d1.10[d1.10$LG==LG,])[1]){
    lines(x=c(0.95,1.05),y=c(max_length-d1.10$Size[d1.10$LG==LG][i],max_length-d1.10$Size[d1.10$LG==LG][i]),col="dark green")
  }
  for(i in 1:dim(d2.15[d2.15$LG==LG,])[1]){
    lines(x=c(1.95,2.05),y=c(max_length-d2.15$flipped_Size[d2.15$LG==LG][i],max_length-d2.15$flipped_Size[d2.15$LG==LG][i]),col="dark blue")
    
    if(d2.15$Marker[d2.15$LG==LG][i] %in% d1.10$Marker[d1.10$LG==LG]){
      lines(x=c(1.15,1.85),y=c(max_length-d1.10$Size[d1.10$LG==LG & as.character(d1.10$Marker)==as.character(d2.15$Marker[d2.15$LG==LG][i])],max_length-d2.15$flipped_Size[d2.15$LG==LG][i]),lty=1,col="black")
    }
  }
}

par(mfrow=c(4,5))
for(i in 1:19){
  plotCompareMap(d1.10,d2.15,i,400)
}


library(gdata)

for (i in 1:19){
  assign(paste("d1.10_LG",i,sep=""),drop.levels(d1.10[d1.10$LG==i,]))
  assign(paste("d2.15_LG",i,sep=""),drop.levels(d2.15[d2.15$LG==i,]))
}


### Flip d2.15 LGs with opostit orientation from d1.10
d2.15_LG1$flipped_Size<-sapply(d2.15_LG1$Size,function(x) {max(d2.15_LG1$Size)-d2.15_LG1$Size})[,1]
d2.15_LG2$flipped_Size<-sapply(d2.15_LG2$Size,function(x) {max(d2.15_LG2$Size)-d2.15_LG2$Size})[,1]
d2.15_LG4$flipped_Size<-sapply(d2.15_LG4$Size,function(x) {max(d2.15_LG4$Size)-d2.15_LG4$Size})[,1]
d2.15_LG5$flipped_Size<-sapply(d2.15_LG5$Size,function(x) {max(d2.15_LG5$Size)-d2.15_LG5$Size})[,1]
d2.15_LG7$flipped_Size<-sapply(d2.15_LG7$Size,function(x) {max(d2.15_LG7$Size)-d2.15_LG7$Size})[,1]
d2.15_LG15$flipped_Size<-sapply(d2.15_LG15$Size,function(x) {max(d2.15_LG15$Size)-d2.15_LG15$Size})[,1]
d2.15_LG17$flipped_Size<-sapply(d2.15_LG17$Size,function(x) {max(d2.15_LG17$Size)-d2.15_LG17$Size})[,1]
d2.15_LG19$flipped_Size<-sapply(d2.15_LG19$Size,function(x) {max(d2.15_LG19$Size)-d2.15_LG19$Size})[,1]

assign("LG1",list(d1.10_LG1[,c(2,3)],d2.15_LG1[,c(2,4)]))
assign("LG2",list(d1.10_LG2[,c(2,3)],d2.15_LG2[,c(2,4)]))
assign("LG3",list(d1.10_LG3[,c(2,3)],d2.15_LG3[,c(2,3)]))
assign("LG4",list(d1.10_LG4[,c(2,3)],d2.15_LG4[,c(2,4)]))
assign("LG5",list(d1.10_LG5[,c(2,3)],d2.15_LG5[,c(2,4)]))
assign("LG6",list(d1.10_LG6[,c(2,3)],d2.15_LG6[,c(2,3)]))
assign("LG7",list(d1.10_LG7[,c(2,3)],d2.15_LG7[,c(2,4)]))
assign("LG8",list(d1.10_LG8[,c(2,3)],d2.15_LG8[,c(2,3)]))
assign("LG9",list(d1.10_LG9[,c(2,3)],d2.15_LG9[,c(2,3)]))
assign("LG10",list(d1.10_LG10[,c(2,3)],d2.15_LG10[,c(2,3)]))
assign("LG11",list(d1.10_LG11[,c(2,3)],d2.15_LG11[,c(2,3)]))
assign("LG12",list(d1.10_LG12[,c(2,3)],d2.15_LG12[,c(2,3)]))
assign("LG13",list(d1.10_LG13[,c(2,3)],d2.15_LG13[,c(2,3)]))
assign("LG14",list(d1.10_LG14[,c(2,3)],d2.15_LG14[,c(2,3)]))
assign("LG15",list(d1.10_LG15[,c(2,3)],d2.15_LG15[,c(2,4)]))
assign("LG16",list(d1.10_LG16[,c(2,3)],d2.15_LG16[,c(2,3)]))
assign("LG17",list(d1.10_LG17[,c(2,3)],d2.15_LG17[,c(2,4)]))
assign("LG18",list(d1.10_LG18[,c(2,3)],d2.15_LG18[,c(2,3)]))
assign("LG19",list(d1.10_LG19[,c(2,3)],d2.15_LG19[,c(2,4)]))

save.image(".RData")