#required packages
library(tidyverse)
install.packages("caret")
install.packages("caret", dependencies=c("Depends", "Suggests"))
library(caret)#
install.packages('Boruta')
library(Boruta)

#input directory for metrics files
metricsDir = "Input/Metrics"
metrics <- data.frame()
met1 <- read.table(list.files(metricsDir,full.names = T)[1],header=T,fill=T)
met1 <- met1 %>% select(Library,Postfiltering_reads_aligned,Sequencing_for_0.05x_cov, Background,Reads_per_Mb,Quality,Sequencing_for_0.05x_cov, Mode_insert_size, Mode_GC, Percent_WC, Naive_coverage, Duplication_rate)
names(met1)[names(met1) == "Postfiltering_reads_aligned"] <- "Reads"
#met1$source=1
met1$Reads <- as.integer(met1$Reads)

met2 <- read.table(list.files(metricsDir,full.names = T)[2],header=T,fill=T)
met2 <- met2 %>% select(Library,Postfiltering_reads_aligned,Sequencing_for_0.05x_cov, Background,Reads_per_Mb,Quality,Sequencing_for_0.05x_cov, Mode_insert_size, Mode_GC, Percent_WC, Naive_coverage, Duplication_rate)
names(met2)[names(met2) == "Postfiltering_reads_aligned"] <- "Reads"
#met2$source=2
met2$Reads <- as.integer(met2$Reads)

met3 <- read.table(list.files(metricsDir,full.names = T)[3],header=T,fill=T)
met3 <- met3 %>% select(Library,Reads_aligned_postfiltering,Sequencing_for_5_percent_coverage, Background,Reads_per_Mb,Quality,Sequencing_for_5_percent_coverage, Mode_insert_size, Mode_GC, Percent_WC, Coverage, Duplication_rate)
names(met3)[names(met3) == "Reads_aligned_postfiltering"] <- "Reads"
names(met3)[names(met3) == "Sequencing_for_5_percent_coverage"] <- "Sequencing_for_0.05x_cov"
names(met3)[names(met3) == "Coverage"] <- "Naive_coverage"
#met3$source=3
met3$Reads <- as.integer(met3$Reads)

met4 <- read.table(list.files(metricsDir,full.names = T)[4],header=T,fill=T)
met4 <- met4 %>% select(Library,Postfiltering_reads_aligned,Background,Reads_per_Mb,Quality, Mode_insert_size, Mode_GC, Percent_WC,Naive_coverage, Duplication_rate)
names(met4)[names(met4) == "Postfiltering_reads_aligned"] <- "Reads"
met4$Sequencing_for_0.05x_cov  = NA
#met4$source=4
met4$Reads <- as.integer(met4$Reads)

met5 <- read.table(list.files(metricsDir,full.names = T)[5],header=T,fill=T)
met5 <- met5 %>% select(Library,Postfiltering_reads_aligned,Sequencing_for_0.05x_cov, Background,Reads_per_Mb,Quality,Sequencing_for_0.05x_cov, Mode_insert_size, Mode_GC, Percent_WC, Naive_coverage, Duplication_rate)
names(met5)[names(met5) == "Postfiltering_reads_aligned"] <- "Reads"
met5$Reads <- as.integer(met5$Reads)
#met5$source=5

#write.table(paste0(metrics$Library,".trimmed.mdup.bam"), "filesInMet.txt",sep="\t",quote=F,row.names = F,col.names = F)
#row binding all metrics files
metrics <- rbind(met1,met2,met3,met5)
metrics$Library=paste0(metrics$Library,".trimmed.mdup.bam")


extraMet <- read.table("../collectLibraryStats.txt",header=T,fill=T)
metrics=merge(metrics,extraMet,by.x="Library",by.y="file")


message("There are a total of ", nrow(metrics)," libraries in all metrics files.")


#Convert percentages to strings and then string splitting, then back to numeric
metrics$Background <- as.character(metrics$Background)
metrics$Background<-strsplit(metrics$Background,"%")
metrics$Background<-as.numeric(metrics$Background)
metrics$Mode_GC <- as.character(metrics$Mode_GC)
metrics$Mode_GC<-strsplit(metrics$Mode_GC,"%")
metrics$Mode_GC<-as.numeric(metrics$Mode_GC)
metrics$Percent_WC <- as.character(metrics$Percent_WC)
metrics$Percent_WC<-strsplit(metrics$Percent_WC,"%")
metrics$Percent_WC<-as.numeric(metrics$Percent_WC)
metrics$Naive_coverage <- as.character(metrics$Naive_coverage)
metrics$Naive_coverage<-strsplit(metrics$Naive_coverage,"%")
metrics$Naive_coverage<-as.numeric(metrics$Naive_coverage)
metrics$Duplication_rate <- as.character(metrics$Duplication_rate)
metrics$Duplication_rate<-strsplit(metrics$Duplication_rate,"%")
metrics$Duplication_rate<-as.numeric(metrics$Duplication_rate)
metrics$Sequencing_for_0.05x_cov <- as.character(metrics$Sequencing_for_0.05x_cov)
metrics$Sequencing_for_0.05x_cov[metrics$Sequencing_for_0.05x_cov==">=2"] <- 2
metrics$Sequencing_for_0.05x_cov[is.na(metrics$Sequencing_for_0.05x_cov)] <- "2"
metrics$Sequencing_for_0.05x_cov <- as.numeric(metrics$Sequencing_for_0.05x_cov )




#convert NAs to 0s
metrics$Reads[is.na(metrics$Reads)] <- 0
metrics$Background[is.na(metrics$Background)] <- 0
metrics$Naive_coverage[is.na(metrics$Naive_coverage)] <- 0
metrics$Duplication_rate[is.na(metrics$Duplication_rate)] <- 0

#annotate empty libraries
for (row in 1:nrow(metrics)){
	#print(metrics[row,]$Reads)
	#row=3
	if (metrics[row,]$Reads<=5000){
		metrics[row,]$Quality<- "shallow"
		#print("empty")
	} #else{print("*")}
}

#metrics=na.exclude(metrics)
metrics=select(metrics,-c(Library))

#ploting
dataset<- metrics
x <- dataset[,c(1:4,6:14)]
y <- dataset[,5]
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot="box",scales =scales)



# create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(metrics$Quality, p=0.80, list=FALSE)
# select 20% of the data for validation
validation <- metrics[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset <- metrics[validation_index,]
sapply(dataset, class)
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
dataset=na.exclude(dataset)

# a) linear algorithms
set.seed(7)
fit.lda <- train(Quality~., data=dataset, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
set.seed(7)
fit.cart <- train(Quality~., data=dataset, method="rpart", metric=metric, trControl=control)
# kNN
set.seed(7)
fit.knn <- train(Quality~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
set.seed(7)
fit.svm <- train(Quality~., data=dataset, method="svmRadial", metric=metric, trControl=control)
# Random Forest
set.seed(7)
fit.rf <- train(Quality~., data=dataset, method="rf", metric=metric, trControl=control)
set.seed(7)
fit.gbm <- train(Quality~., data=dataset, method="gbm", metric=metric, trControl=control)
set.seed(7)
fit.adaboost <- train(Quality~., data=dataset, method="adaboost", metric=metric, trControl=control)


results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf,gbm=fit.gbm))
summary(results)
dotplot(results)
print(fit.rf)

#importance
fit.rf.Imp <- varImp(fit.rf, scale = FALSE)
plot(fit.rf.Imp)
boruta_output <- Boruta(Quality ~ ., data=na.omit(dataset), doTrace=0)
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
head(imps2[order(-imps2$meanImp), ])  # descending sort

validation=na.omit(validation)
predictions <- predict(fit.rf, na.omit(validation))
confusionMatrix(predictions, validation$Quality)
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")
str(metrics$Background)

default_knn_mod = train(
	Quality ~ .,
	data = dataset,
	method = "knn",
	trControl = trainControl(method = "cv", number = 5),
	preProcess = c("center", "scale"),
	tuneGrid = expand.grid(k = seq(1, 101, by = 2)))
plot(default_knn_mod)
plot(fit.rf)

featurePlot(x = metrics[, c("complexity","spikiness","entropy","Background","Reads")],
			y = metrics$Quality,
			plot = "pairs",
			auto.key = list())

pairs(metrics[, c("complexity","spikiness","entropy","Background","Reads")],
	  col=metrics$Quality)
pca(metrics[2:3])
BiocManager::install("M3C")
library(M3C)


