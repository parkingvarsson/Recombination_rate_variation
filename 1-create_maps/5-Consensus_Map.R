load("/mnt/picea/projects/u2016004/Genetic-map/Round03/.RData")

library(LPmerge)
args <- commandArgs(trailingOnly = TRUE)
LG<-get(paste("LG",args[1],sep=""))

Consensus <- LPmerge(LG,max.interval=1:10)

out_name<-paste("LG",args[1],"_Consensus",sep="")
assign(out_name,Consensus)
save(list=out_name,file=paste("LG",args[1],"_Consensus.RData",sep=""))