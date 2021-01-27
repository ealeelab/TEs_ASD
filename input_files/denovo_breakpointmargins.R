callset<-read.table("p1_HighConfInsertions_calls.txt",sep="\t",header=FALSE)

#if onesided then use pos and add -40 +40 margin, but otherwise then take midpoint and add margin

margin=50

callset$s <- ifelse((callset$V3==-1 | callset$V4==-1), callset$V2-margin, floor((callset$V3+callset$V4)/2)-margin)
callset$e <- ifelse((callset$V3==-1 | callset$V4==-1), callset$V2+margin, floor((callset$V3+callset$V4)/2)+margin)  

#if negative make pos zero start and end margin
callset$s[callset$s<0]<-1
callset$e[callset$e<margin]<-margin

callsetfinal<-data.frame(callset$V1,callset$s,callset$e)

write.table(callsetfinal,"p1_HighConfInsertions_calls.bed", sep="\t", quote=FALSE,col.names = FALSE, row.names = FALSE)
#  save new coordinates with original coordinates
write.table(callset,"p1_HighConfInsertions_calls_Original.txt", sep="\t", quote=FALSE,col.names = FALSE, row.names = FALSE)

quit(save="no")


