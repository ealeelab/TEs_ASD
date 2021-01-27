clipped_TE<-read.table("p1_HighConfInsertions_noDISCparent_noTECLIPparent.bed",sep="\t",header=FALSE)
clipped_all<-read.table("p1_HighConfInsertions_noDISCparent_noALLCLIPparent.bed",sep="\t",header=FALSE)

#input example: 10	24858127	24858227

clipped_merged<-merge(clipped_TE,clipped_all,by=c("V1","V2","V3"))
clipped_merged<-unique(clipped_merged)

write.table(clipped_merged,"p1_HighConfInsertions_noDISCparent_noCLIPparent.bed", sep="\t", quote=FALSE,col.names = FALSE, row.names = FALSE)
quit(save="no")

