##### Figure 2 #####


library("data.table")
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("cowplot")

df<-fread("~/Dropbox/Linked_selection/recomb_vs/All_features_recomb_meth.txt",head=T)[,c(1:15)]
div<-fread("/Users/pron0005/Dropbox/Linked_selection/recomb_vs/substitutions/SRR1569781.HC-p2.chr.sorted.fixed_substitutions.neutral.all.1Mbp_window.snpden")
df<-merge(df,div,by=c("chromosome","win_start"))
df$chr<-as.numeric(substr(df$chromosome,4,6))
df<-setorder(df,chr,win_start)
df$win_start<-df$win_start/1e6

col1<-"#4b7c6c"
col2<-"#67b49b"
col3<-"#91ccb8"
col4<-"#b4acac"
col5<-"#fb8c62"
col6<-"#e39c81"
col7<-"#e8c2b1"

df1<-subset(df,chr==1)
df2<-subset(df,chr==5)
p1a <- ggplot(df1,aes(x=win_start,y=gen_slidingwindow)) + geom_area(fill=col1) + xlab("") + coord_cartesian(ylim=c(0,25)) + ylab("cM/Mb") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2a <- ggplot(df2,aes(x=win_start,y=gen_slidingwindow)) + geom_area(fill=col1) + xlab("") + coord_cartesian(ylim=c(0,25)) + ylab("cM/Mb") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p1b <- ggplot(df1,aes(x=win_start,y=seq_slidingwindow)) + geom_area(fill=col2) + coord_cartesian(ylim=c(0,25)) + xlab("") + ylab("cM/Mb") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2b <- ggplot(df2,aes(x=win_start,y=seq_slidingwindow)) + geom_area(fill=col2) + coord_cartesian(ylim=c(0,25)) + xlab("") + ylab("cM/Mb") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p1c <- ggplot(df1,aes(x=win_start,y=neut_div)) + geom_area(fill=col3) + xlab("") + ylab("Nucl. div") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.00025,by=0.0001), labels = scales::scientific)
p2c <- ggplot(df2,aes(x=win_start,y=neut_div)) + geom_area(fill=col3) +  xlab("") + ylab("Nucl. div") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.00025,by=0.0001), labels = scales::scientific)

p1d <- ggplot(df1,aes(x=win_start,y=total_neut_subs_density)) + geom_area(fill=col4) +  xlab("") + ylab("Divergence") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,6,by=2))
p2d <- ggplot(df2,aes(x=win_start,y=total_neut_subs_density)) + geom_area(fill=col4) +  xlab("") + ylab("Divergence") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,6,by=2))

p1e <- ggplot(df1,aes(x=win_start,y=gene_cov)) + geom_area(fill=col5) + xlab("")  + ylab("Gene density") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.5,by=0.25))
p2e <- ggplot(df2,aes(x=win_start,y=gene_cov)) + geom_area(fill=col5) + xlab("")  + ylab("Gene density") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.5,by=0.25))

p1f <- ggplot(df1,aes(x=win_start,y=repeat_cov)) + geom_area(fill=col6) +  xlab("") + ylab("Repeat density") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.25,by=0.1))
p2f <- ggplot(df2,aes(x=win_start,y=repeat_cov)) + geom_area(fill=col6) +  xlab("") + ylab("Repeat density") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks = seq(0,0.25,by=0.1))

p1g <- ggplot(df1,aes(x=win_start,y=CpG_mean)) + geom_area(fill=col7) + coord_cartesian(ylim=c(0,60)) +  xlab("Chr 1 (Mb)") + ylab("CpG meth %") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2g <- ggplot(df2,aes(x=win_start,y=CpG_mean)) + geom_area(fill=col7) + coord_cartesian(ylim=c(0,60)) +  xlab("Chr 5 (Mb)") + ylab("CpG meth %") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


plot_grid(p1a,p2a,p1b,p2b,p1c,p2c,p1d,p2d,p1e,p2e,p1f,p2f,p1g,p2g,ncol=2,align="v",rel_widths = c(1.92,1))
dev.off()

p1 <- ggplot(df2,aes(x=win_start)) + geom_line(aes(y=gen_slidingwindow)) + xlab("") + ylab("cM/Mb") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
