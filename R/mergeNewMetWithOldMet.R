met <- read.table("../../Downloads/metrics_summary-2.txt",header=T,fill=T)
met <- met %>% select(Library,Postfiltering_reads_aligned, Background,Reads_per_Mb,Sequencing_for_0.05x_cov, Mode_insert_size, Mode_GC, Percent_WC, Naive_coverage, Duplication_rate)
met$Library=paste0(met$Library,".trimmed.mdup.bam")

#Convert percentages to strings and then string splitting, then back to numeric
met$Background <- as.character(met$Background)
met$Background<-strsplit(met$Background,"%")
met$Background<-as.numeric(met$Background)
met$Mode_GC <- as.character(met$Mode_GC)
met$Mode_GC<-strsplit(met$Mode_GC,"%")
met$Mode_GC<-as.numeric(met$Mode_GC)
met$Percent_WC <- as.character(met$Percent_WC)
met$Percent_WC<-strsplit(met$Percent_WC,"%")
met$Percent_WC<-as.numeric(met$Percent_WC)
met$Naive_coverage <- as.character(met$Naive_coverage)
met$Naive_coverage<-strsplit(met$Naive_coverage,"%")
met$Naive_coverage<-as.numeric(met$Naive_coverage)
met$Duplication_rate <- as.character(met$Duplication_rate)
met$Duplication_rate<-strsplit(met$Duplication_rate,"%")
met$Duplication_rate<-as.numeric(met$Duplication_rate)
met$Sequencing_for_0.05x_cov <- as.character(met$Sequencing_for_0.05x_cov)
met$Sequencing_for_0.05x_cov[met$Sequencing_for_0.05x_cov==">=2"] <- 2
met$Sequencing_for_0.05x_cov[is.na(met$Sequencing_for_0.05x_cov)] <- "2"
met$Sequencing_for_0.05x_cov <- as.numeric(met$Sequencing_for_0.05x_cov )


#convert NAs to 0s
met$Reads[is.na(met$Reads)] <- 0
met$Background[is.na(met$Background)] <- 0
met$Naive_coverage[is.na(met$Naive_coverage)] <- 0
met$Duplication_rate[is.na(met$Duplication_rate)] <- 0

#annotate empty libraries
for (row in 1:nrow(met)){
	#print(met[row,]$Reads)
	#row=3
	if (met[row,]$Reads<=5000){
		met[row,]$Quality<- "shallow"
		#print("empty")
	} #else{print("*")}
}

extraMet <- read.table("../JAN_2021_BAM.txt",header=T,fill=T)
metrics=merge(met,extraMet,by.x="Library",by.y="file")


metrics$predictions <- predict(fit.rf, metrics)

goodPred <- filter(metrics,predictions=="good")
nrow(goodPred)

write.table(goodPred$Library,"goodPred.txt",sep="\t",quote=F,col.names = F,row.names = F)



