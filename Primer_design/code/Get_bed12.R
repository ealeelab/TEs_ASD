
# Code to convert calls into bed12
# Change directory and Input file paths

setwd("")
callset<-read.table("",sep="\t",header=TRUE)

# Uses the following input format:
# FAMILY	SAMPLE_ID	TYPE	CHR	START	END	PHASE	FILTER_DENOVO	CONFIDENCE	Gene	Pli	SAFARI
#HG002	NA24385	L1	chr14	30681602	30681617	NA	two_side_tprt_both	2	NA	NA	NA
#HG002	NA24385	L1	chr5	90154952	90154962	NA	two_side_tprt_both	2	NA	NA	NA
#HG002	NA24385	L1	chr3	38584574	38584591	NA	two_side_tprt_both	2	NA	NA	NA


# coordinates only

#if onesided then use pos and add -/+ margin, but otherwise then take midpoint and add -/+ margin

margin=1

# If no breakpoint, do -/+ existing breakpoint

callset$START <- ifelse((callset$START==-1 ), callset$END-margin , callset$START)
callset$END <- ifelse((callset$END==-1), callset$START+margin , callset$END)  

# make sure start is smaller
callset$bedstart <- ifelse((callset$START < callset$END ), callset$START , callset$END)
callset$bedend <- ifelse((callset$END > callset$START ), callset$END , callset$START)  

# get bed12 

# we want
# chr	start-800	end+800	FAMILY+SAMPLE_ID_TYPE	0	+	start	start	NA	2 	700,700	0,end-700

callset$bedstart <- callset$bedstart - 800
callset$bedend  <- callset$bedend +800

callset$bedblockStarts_endpoint<-callset$bedend-700
callset$bedblockdistendpoint<-abs(callset$bedstart-callset$bedblockStarts_endpoint)

callset$bedID <- paste(callset$FAMILY,callset$SAMPLE_ID,callset$TYPE,callset$CHR,callset$START,sep = "_")

callset$bedScore <- 0
callset$bedStrand <- "+"
callset$bedthickStart <- callset$START
callset$bedthickEnd <- callset$START
callset$bedRgb <- "NA"
callset$bedblockCount <- 2
callset$bedblockLength <- "700,700"
callset$bedblockStarts <- paste("0",callset$bedblockdistendpoint,sep = ",")

callsetfinal<-data.frame(callset$CHR,callset$bedstart,callset$bedend,callset$bedID,callset$bedScore,callset$bedStrand,callset$bedthickStart,callset$bedthickEnd,callset$bedRgb,callset$bedblockCount,callset$bedblockLength,callset$bedblockStarts)


write.table(callsetfinal,"Insertion_Candidates.bed12", sep="\t", quote=FALSE,col.names = FALSE, row.names = FALSE)

quit(save="no")


