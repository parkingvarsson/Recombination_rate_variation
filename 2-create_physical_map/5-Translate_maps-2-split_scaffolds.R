read.csv("Female_with_bins.csv",head=T)->Female_orig
read.csv("Male_with_bins.csv",head=T)->Male_orig
read.table("v1.1_Potra01-genome.fa.masked.agp",head=F)->scaffold_translate

Female_translated_map<-as.data.frame(matrix(NA,nrow = dim(Female_orig)[1],ncol=dim(Female_orig)[2]))
names(Female_translated_map)<-c("Scaffold ID","scaffold position","LG","genetic position")

for(i in 1:dim(Female_orig)[1]){
  scaffold<-scaffold_translate[scaffold_translate$V6==as.character(Female_orig$Scaffold.ID[i]),]
  for(j in 1:dim(scaffold)[1]){
    if(Female_orig$scaffold.position[i] > scaffold$V7[j] & Female_orig$scaffold.position[i] < scaffold$V8[j]){
      Female_translated_map$`Scaffold ID`[i]<-as.character(scaffold$V1[j])
      Female_translated_map$`scaffold position`[i]<-Female_orig$scaffold.position[i] - scaffold$V7[j] + scaffold$V2[j]
      Female_translated_map$LG[i]<-Female_orig$Chr[i]
      Female_translated_map$`genetic position`[i]<-Female_orig$genetic.position[i]
    }
  }
}

Male_translated_map<-as.data.frame(matrix(NA,nrow = dim(Male_orig)[1],ncol=dim(Male_orig)[2]))
names(Male_translated_map)<-c("Scaffold ID","scaffold position","LG","genetic position")

for(i in 1:dim(Male_orig)[1]){
  scaffold<-scaffold_translate[scaffold_translate$V6==as.character(Male_orig$Scaffold.ID[i]),]
  for(j in 1:dim(scaffold)[1]){
    if(Male_orig$scaffold.position[i] > scaffold$V7[j] & Male_orig$scaffold.position[i] < scaffold$V8[j]){
      Male_translated_map$`Scaffold ID`[i]<-as.character(scaffold$V1[j])
      Male_translated_map$`scaffold position`[i]<-Male_orig$scaffold.position[i] - scaffold$V7[j] + scaffold$V2[j]
      Male_translated_map$LG[i]<-Male_orig$Chr[i]
      Male_translated_map$`genetic position`[i]<-Male_orig$genetic.position[i]
    }
  }
}

write.csv(Male_translated_map,"F1_Male.csv",quote = F,row.names = F)
write.csv(Female_translated_map,"F1_Female.csv",quote = F,row.names = F)