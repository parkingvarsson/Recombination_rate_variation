library(data.table)
#' Read only the F1
F1 <- fread("F1.GT.FORMAT.txt",header=TRUE)
#' Extract the parents
Parents <- F1[,.(CHROM,POS,QUAL,DIST,MAL,FEM,
                 UME_303601_P01_WA05,UME_303601_P01_WA08)]
agree <-Parents[UME_303601_P01_WA08 != "./."][UME_303601_P01_WA05 != "./."][,
  .("fagree"=sum(FEM==UME_303601_P01_WA05),
  "fdisagree"=sum(FEM!=UME_303601_P01_WA05),
  "magree"=sum(MAL==UME_303601_P01_WA08),
  "mdisagree"=sum(MAL!=UME_303601_P01_WA08))]
#' Remove the parents from the original table
F1[,c("UME_303601_P01_WA05","UME_303601_P01_WA08","MAL","FEM") := NULL]

#' Generate a marker index as "Chromosome:Position"
Parents[,marker := .(paste(CHROM, POS, sep=":"))]
F1[,marker := .(paste(CHROM, POS, sep=":"))]
#' Use the generated marker as a table index
setkey(F1, "marker")
setkey(Parents,  "marker")
#' Count the genotypes/markers, a="0/0", b="1/1", ab="0/1", mis¬ßs="./."
#' C++ program for memory efficient, fast counting
F1_segregation <-
  fread(paste("cut -f 1,2,9-772 F1.GT.FORMAT.txt",
              "| ./segregation - extended_sex_rapid_combined_probes.bed"),
        header = TRUE)
setkey(F1_segregation, "marker")

#' Two data.tables with the same key column can be joined as `A[B]`
F1_segregation <-
  F1_segregation[Parents[,.(UME_303601_P01_WA05, UME_303601_P01_WA08,
                            MAL, FEM, marker)]]

F1_segregation <- F1_segregation[nother < 10,]
F1 <- F1[F1_segregation[,marker]]
Parents <- Parents[F1_segregation[,marker]]

##recount wronly assigned genotypes (according to parental genotypes)
##to missing data and set crosstype
df<-F1_segregation
df[,crosstype := character(nrow(df))]
if(sum(df[FEM == "0/0" & MAL == "0/1", b]) > 0)
  warn("You have impossible genotypes. Re-run vcf2gt")
df[FEM == "0/0" & MAL == "0/1", crosstype := "D2.15"]

if(sum(df[FEM == "1/1" & MAL == "0/1", a]) > 0)
  warn("You have impossible genotypes. Re-run vcf2gt")
df[FEM == "1/1" & MAL == "0/1", crosstype := "D2.15"]

if(sum(df[FEM == "0/1" & MAL == "0/0", b]) > 0)
  warn("You have impossible genotypes. Re-run vcf2gt")
df[FEM == "0/1" & MAL == "0/0", crosstype := "D1.10"]

if(sum(df[FEM == "0/1" & MAL == "1/1", a]) > 0)
  warn("You have impossible genotypes. Re-run vcf2gt")
df[FEM == "0/1" & MAL == "1/1", crosstype := "D1.10"]

df[FEM == "0/1" & MAL == "0/1", crosstype := "B3.7"]

##discard all markers with >20% missing data (0.2*(ncol(F1_final)-3))
df <- df[miss < 0.2 * (ncol(F1) - 5),]
##discard uninformative parents and those outside probe regions
df <- df[MAL != "./." & FEM != "./."]
df <- df[!is.na(probe)]

##calculate segregation distortion
double_het_markers <- df[FEM == "0/1" & MAL == "0/1",]
hom_ref_markers    <- df[FEM == "0/0" | MAL == "0/0",]
hom_alt_markers    <- df[FEM == "1/1" | MAL == "1/1",]

double_het_markers[,pvalue := chisq.test(x=c(a, b, ab),
                                         p=c(0.25,0.25,0.5))[[3]], by = marker]
hom_ref_markers[,pvalue := chisq.test(c(a, ab, b))[[3]], by = marker]
hom_alt_markers[,pvalue := chisq.test(c(ab, b))[[3]], by = marker]

df[double_het_markers[,marker], pvalue := double_het_markers[,pvalue]]
df[hom_ref_markers[,marker], pvalue := hom_ref_markers[,pvalue]]
df[hom_alt_markers[,marker], pvalue := hom_alt_markers[,pvalue]]

## create a file with informative markers (pvalue >0.005)
informative_markers<-df[F1[,.(QUAL,DIST,marker)]][!is.na(probename)][pvalue >= 0.005]
informative_markers[,order := order(probename,QUAL,decreasing = TRUE)]
setkey(informative_markers,order)
write.table(informative_markers,file = "F1-informative_markers-all.txt",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
informative_markers<-informative_markers[!duplicated(probename)]
informative_markers<-informative_markers[,.(CHROM,POS,probename)]
write.table(informative_markers,"F1-informative_markers.txt",
            quote=FALSE,row.names=FALSE,sep="\t" )

##To_Onemap
info_markers_F1 <- F1[df[pvalue >= 0.005][!duplicated(probename),
                         .(marker,crosstype)]]

#F1_segregation_final <- merge(F1_segregation,info_markers_F1,by="marker")

#info_markers_F1[,1:3,with = FALSE]

cols <- grep("marker|UME", names(info_markers_F1))

info_markers_F1 <- info_markers_F1[crosstype != ""]

F1_Onemap <-
  data.table(marker = info_markers_F1[,paste("*", marker, sep ="")],
             crosstype = info_markers_F1[,crosstype],
             GT = info_markers_F1[,cols,with = FALSE][,.(GT = paste(.SD,collapse=",")),
                                                      by = marker][,GT])

transform <- function(crosstype, sequence)
{
  v <- unlist(strsplit(sequence, ","))
  if(crosstype == "B3.7"){
    v <- gsub("1/1","b",v)
  }
  else {
    v <- gsub("1/1","a",v)
  }
  v <- gsub("0/0","a",v)
  v <- gsub("0/1","ab",v)
  v <- gsub("./.","-",v)
  return(paste(v,collapse=","))
}

F1_Onemap <- F1_Onemap[,.(crosstype, GT = transform(crosstype, GT)),by = marker]

file.out<-file("F1_Onemap.txt","w")
cat(c(ncol(F1)-5,dim(F1_Onemap)[1],0),sep="\t",file=file.out)
cat("\n",file=file.out)
write.table(F1_Onemap,col.names=F, row.names=F, quote=F,sep="\t",file=file.out,append=T)
flush(file.out)
close(file.out)
