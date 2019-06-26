###Figure 3###
####Smoothed/Windowed####
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")
### Comparison ###
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")

## Create a correlation matrix ##
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")
neutral_substitutions <- read.table("/path/to/substitutions/SRR1569781.HC-p2.chr.sorted.fixed_substitutions.neutral.all.1Mbp_window.snpden", header=T)
All_features_recomb$old_subs_density <- neutral_substitutions$old_neut_subs_density
All_features_recomb$new_subs_density <- neutral_substitutions$new_neut_subs_density
All_features_recomb <- All_features_recomb[,c(1,2,3,4,5,6,9,8,18,17,7,10,12,14,16,11,13,15)]
All_features_recomb$CpG_mean <- All_features_recomb$CpG_mean/100
All_features_recomb$CHG_mean <- All_features_recomb$CHG_mean/100
All_features_recomb$CHH_mean <- All_features_recomb$CHH_mean/100
names(All_features_recomb) <- c("chromosome", "win_start", "win_end", "Gen. based", "Seq. based", "Gene dens.", "GC-content", "Neutral div.", "New sub. dens.", "Old sub. dens.", "Repeat dens.", "CpG", "CHG", "CHH")
correlation_matrix <- round(cor(All_features_recomb[,c(4:14)], use="pairwise", method = "spearman"), 3)

library(pheatmap)
m<-correlation_matrix
labels<-c("LMB estimates", "LDB estimates", "Gene dens.", "GC-content", "Neutral div.", "New subs. dens.", "Old subs. dens.", "Repeat dens.", "CpG", "CHG", "CHH")
paletteLength <- 100
mapcolors <- colorRampPalette(c("blue","white","red"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(m)/paletteLength, max(m), length.out=floor(paletteLength/2)))
pheatmap(m,cluster_rows=FALSE,cluster_cols=FALSE, breaks = myBreaks, color = mapcolors,
         main = "",#gaps_row = c(5,8),gaps_col = c(5,8),
         display_numbers = T, #matrix(ifelse(is.na(m)|m<1,"",m), nrow(m)),
         na_col = "white",labels_col = labels,
         labels_row = labels,fontsize_number = 11, fontsize = 10)
dev.off()

###Figure S8###
####Smoothed/Windowed####
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")
### Comparison ###
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")

## Make chromosome column into numeric ##
temp <- gsub("chr([0-9]+)", "\\1", All_features_recomb$chromosome)
All_features_recomb$chromosome <- as.integer(temp)
rm(temp)

## Draw the plot ##
# Call the "manhattan_plot_modified.R" before running this #
par(mfrow=c(11,1), mar=c(1,5,0.1,0.1), oma=c(3,0.1,0.1,0.1))

manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="gen_slidingwindow", logp=F, ylabel="Gen. rec.",
          type="l", lwd=2, ymin=0, ymax=(max(All_features_recomb$gen_slidingwindow, na.rm=T)+(max(All_features_recomb$gen_slidingwindow, na.rm=T)/25)))
title(ylab="GMB est.", line=4)

manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="seq_slidingwindow", logp=F, ylabel="Seq. rec.",
          type="l", lwd=2, ymin=0, ymax=(max(All_features_recomb$gen_slidingwindow, na.rm=T)+(max(All_features_recomb$gen_slidingwindow, na.rm=T)/25)))
title(ylab="LDB est.", line=4)

manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="gene_cov", logp=F, ylabel="Gene cov.",
          type="l", lwd=2, ymin=0, ymax=0.6) #(max(All_features_recomb$gene_cov, na.rm=T)+(max(All_features_recomb$gene_cov, na.rm=T)/25)))
title(ylab="Gene dens.", line=4)

manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="gccontent", logp=F, ylabel="GC-content",  type="l", lwd=2, ymin=0.31, ymax=0.37) #(max(All_features_recomb$gccontent, na.rm=T)+(max(All_features_recomb$gccontent, na.rm=T)/25)))
title(ylab="GC-cont.", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="neut_div", logp=F, ylabel="Neut. div.",  type="l", lwd=2, ymin=0, ymax=(max(All_features_recomb$neut_div, na.rm=T)+(max(All_features_recomb$neut_div, na.rm=T)/25)))
title(ylab="Neut. div.", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="repeat_cov", logp=F, ylabel="Repeat cov.",
          type="l", lwd=2, ymin=0, ymax=0.3) #(max(All_features_recomb$repeat_cov, na.rm=T)+(max(All_features_recomb$repeat_cov, na.rm=T)/25)))
title(ylab="Repeat dens.", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="CpG_mean", logp=F, ylabel="Methylation %",  type="l", lwd=2, ymin=0, ymax=60) #=(max(All_features_recomb$meth_mean, na.rm=T)+(max(All_features_recomb$meth_mean, na.rm=T)/25)))
title(ylab="CpG %", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="CHG_mean", logp=F, ylabel="Methylation %",  type="l", lwd=2, ymin=0, ymax=40) #=(max(All_features_recomb$meth_mean, na.rm=T)+(max(All_features_recomb$meth_mean, na.rm=T)/25)))
title(ylab="CHG %", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="CHH_mean", logp=F, ylabel="Methylation %",  type="l", lwd=2, ymin=0, ymax=10) #=(max(All_features_recomb$meth_mean, na.rm=T)+(max(All_features_recomb$meth_mean, na.rm=T)/25)))
title(ylab="CHH %", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="old_subs_density", logp=F, ylabel="Old subs/kb",  type="l", lwd=2, ymin=0, ymax=22) #=(max(All_features_recomb$meth_mean, na.rm=T)+(max(All_features_recomb$meth_mean, na.rm=T)/25)))
title(ylab="Old s./kb", line=4)


manhattan(All_features_recomb, chr="chromosome", bp="win_start", p="new_subs_density", logp=F, ylabel="New subs/kb",  type="l", createx=T, lwd=2, ymin=0, ymax=2.5) #=(max(All_features_recomb$meth_mean, na.rm=T)+(max(All_features_recomb$meth_mean, na.rm=T)/25)))
title(ylab="New s./kb", line=4)


mtext("Chromosome", side=1, line=2)

dev.off()




###Figure S9###
####Smoothed/Windowed####
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")
### Comparison ###
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")

plot(All_features_recomb$gen_slidingwindow, All_features_recomb$seq_slidingwindow, xlab="gen. recomb.", ylab="seq. recomb.", main="All chromosomes", xlim=c(0,30), ylim=c(0,40), pch=16, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.4), las=1)
library(Hmisc)
#subplot(plot(All_features_recomb$gen_slidingwindow, All_features_recomb$seq_slidingwindow, xlab="", ylab="", pch=16, cex=0.3, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.4), cex.axis=0.5, las=1), c(0.5,5.5), c(25,40), abline(1,1))
inner_plot <- subplot(plot(All_features_recomb$gen_slidingwindow, All_features_recomb$seq_slidingwindow, xlab="", xlim=c(0,30), xaxt="n", ylab="", yaxt="n", pch=16, cex=0.3, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.4), tck=-0.05, las=1), c(0.5,10.5), c(30,40))
op<-par(no.readonly=TRUE)
par(inner_plot)
clip(-100,30,-100,40)
abline(v = 30, col="red")
abline(h = 40, col="red")
axis(2, las=1, tck=-0.05, mgp=c(3, .4, 0), cex.axis=0.5)
axis(1, tck=-0.05, mgp=c(3, .01, 0), cex.axis=0.5)
par(op)
dev.off()

###Figure S7###
boxplot(All_features_recomb$gen_slidingwindow, All_features_recomb$seq_slidingwindow,names = c("Genetic linkage map based", "Sequence based"), col = c("purple", "orange"), pch=19, cex=0.5, outcol="gray", ylim=c(0,50), ylab="Recombination rate (cM/Mb)", las=1)
points(c(mean(All_features_recomb$gen_slidingwindow, na.rm=T),mean(All_features_recomb$gen_slidingwindow, na.rm=T)), col="red", pch=19)
legend("topleft","(x,y)", c("Mean","Outlier"),col=c("red", "gray"), pch=c(19,19), pt.cex=c(1,0.5))
dev.off()

###Table S2###
## Summary statistics ##
All_features_recomb <- read.csv("/path/to/All_features_recomb_meth.txt", sep = "\t")
summary <- summary(All_features_recomb)
