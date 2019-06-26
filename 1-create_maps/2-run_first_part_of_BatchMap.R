library(onemap)

read_outcross_cpp("/mnt/picea/projects/aspseq/u2016004/Genetic-map/F1_Onemap.txt")->F1.out
save.image(".RData")
twopts_F1 <- rf.2pts(F1.out, LOD=8, max.rf=0.35)

mark.all_F1 <- make.seq(twopts_F1, "all")
save.image(".RData")

class(F1.out)[2]<- "outcross"
class(twopts_F1)[2] <- "outcross"

LGs_F1<- group(mark.all_F1, LOD = 12) ##group(input.seq, LOD=NULL, max.rf=NULL)
print(LGs_F1,detailed=F) ## MAKE SURE THAT THE DATA GROUPS INTO 19 LGs, OTHERWISE CHANGE LOD SCORE within group() command

## make LGs into sequences
for (i in 1:LGs_F1$n.groups){
  print(i)
  assign(paste("LG",i,sep="_"),make.seq(LGs_F1,i))
}

for(f in 1:LGs_F1$n.groups)
{
  lg <- get(paste("LG_",f,sep=""))
  outa <- paste("LG_",f,"_p1",sep="")
  outb <- paste("LG_",f,"_p2",sep="")
  assign(outa, make.seq(twopts_F1,lg$seq.num[F1.out$segr.type[lg$seq.num] != "D2.15"]))
  assign(outb, make.seq(twopts_F1,lg$seq.num[F1.out$segr.type[lg$seq.num] != "D1.10"]))
}

save.image(".RData")

for(f in 1:LGs_F1$n.groups)
{
  lg <- get(paste("LG_",f,sep=""))
  outa <- paste("LG_",f,"_p1",sep="")
  outb <- paste("LG_",f,"_p2",sep="")
  cat(paste(outa, length(get(outa)$seq.num),"\n"))
  cat(paste(outb, length(get(outb)$seq.num),"\n"))
}
