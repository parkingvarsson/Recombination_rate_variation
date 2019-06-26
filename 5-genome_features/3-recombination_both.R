##### Script for merging the two recombination rates together into one file #####

consensus_gen <- read.table("/path/to/MareyMap/outputs/Consensus.MM.output.txt", header=TRUE)
consensus_seq <- read.table("/path/to/MareyMap/outputs/sequence.point.MM_output.txt", header=TRUE)

### Merge ###

consensus_both<-merge(consensus_gen[,c(2,4,5,7)],consensus_seq[,c(2,4,5,7)],by=c("map","phys"),all=T)
names(consensus_both)<-c("chromosome","position","gen_gen","gen_slidingwindow","seq_gen","seq_slidingwindow")

write.table(consensus_both, "/path/to/wherever/you/want/both_maps.slidingwindow.txt", sep = "\t", quote = FALSE, row.names = FALSE)
