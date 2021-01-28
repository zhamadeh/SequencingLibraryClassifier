################### Required packages ##################
install.packages("caret", dependencies=c("Depends", "Suggests"))
install.packages('Boruta')
install.packages('RANN')
library(tidyverse)
library(caret)
library(Boruta)
library(tidyr)
library(RANN)

################### Assembly of training set ##################

#consists of two data sources: standard library metrics + custom new metrics from `collectLibraryMetrics.R`
metricsDir = "Input/Metrics"
metrics <- data.frame()

for (file in list.files(metricsDir,full.names = T)){
	met <- read.table(file,header = T,fill=T,row.names = NULL)
	#message(basename(file), ' has ',nrow(met)," liibraries which is added to ",nrow(metrics), " making ")
	if ("Library" %in% colnames(met) && "Postfiltering_reads_aligned" %in% colnames(met) && "Sequencing_for_0.05x_cov" %in% colnames(met) && "Background" %in% colnames(met) && "Reads_per_Mb" %in% colnames(met) && "Quality" %in% colnames(met) && "Sequencing_for_0.05x_cov" %in% colnames(met) && "Mode_insert_size" %in% colnames(met) && "Mode_GC" %in% colnames(met) && "Percent_WC" %in% colnames(met) && "Naive_coverage" %in% colnames(met) && "Duplication_rate" %in% colnames(met)){
		met <- met %>% select(Library,Postfiltering_reads_aligned,Sequencing_for_0.05x_cov, Background,Reads_per_Mb,Quality,Sequencing_for_0.05x_cov, Mode_insert_size, Mode_GC, Percent_WC, Naive_coverage, Duplication_rate)
		names(met)[names(met) == "Postfiltering_reads_aligned"] <- "Reads"
	} else if ("Reads_aligned_postfiltering" %in% colnames(met) && "Sequencing_for_5_percent_coverage" %in% colnames(met)){
		met <- met %>% select(Library,Reads_aligned_postfiltering,Sequencing_for_5_percent_coverage, Background,Reads_per_Mb,Quality, Mode_insert_size, Mode_GC, Percent_WC, Coverage, Duplication_rate)
		names(met)[names(met) == "Reads_aligned_postfiltering"] <- "Reads"
		names(met)[names(met) == "Sequencing_for_5_percent_coverage"] <- "Sequencing_for_0.05x_cov"
		names(met)[names(met) == "Coverage"] <- "Naive_coverage"
	}
	met$Reads <- as.integer(met$Reads)
	met$Library=paste0(met$Library,".trimmed.mdup.bam")
	metrics <- rbind(met,metrics)
	#message(nrow(metrics)," liibraries")
}
message("There are ",length(metrics$Library)," libraries in the base metrics file with ",ncol(metrics)," features")

#Addition of new metrics
extraMetricsDir = "Input/NewMetrics/"
extrametrics <- data.frame()

for (file in list.files(extraMetricsDir,full.names = T)){
	met <- read.table(file,header = T,fill=T,row.names = NULL)
	#message(basename(file), ' has ',nrow(met)," liibraries which is added to ",nrow(extrametrics), " making ")
	extrametrics <- rbind(met,extrametrics)
	#message(nrow(extrametrics)," liibraries")
}

message("There are ",length(extrametrics$file)," libraries in the new metrics file with ",ncol(extrametrics)," features")
metrics=merge(metrics,extrametrics,by.x="Library",by.y="file")
message("All together, there are now ",nrow(metrics), " libraries with ", ncol(metrics), " features")

#Convert percentages to strings and then string splitting, then back to numeric
for (colname in colnames(metrics)){
	if (grepl( "%", as.character(metrics[,colname][50]), fixed = TRUE)){
		metrics[,colname] <- as.character(metrics[,colname])
		metrics[,colname] <- str_remove(metrics[,colname],"%")
		metrics[,colname]<-as.numeric(metrics[,colname])
	} else if (colname=="Sequencing_for_0.05x_cov"){
		metrics$Sequencing_for_0.05x_cov <- as.character(metrics$Sequencing_for_0.05x_cov)
		metrics$Sequencing_for_0.05x_cov <- str_remove(metrics$Sequencing_for_0.05x_cov,">=")
		metrics$Sequencing_for_0.05x_cov <- as.numeric(metrics$Sequencing_for_0.05x_cov )
		metrics$Sequencing_for_0.05x_cov[is.na(metrics$Sequencing_for_0.05x_cov)] <- 2 #replace_na(metrics$Sequencing_for_0.05x_cov) <- 2
	}
	if (colname == "Percent_WC"|colname == "complexity" ){
		metrics[,colname][is.na(metrics[,colname])] <- 0
	}
}
#annotate empty libraries
for (row in 1:nrow(metrics)){
	if (metrics[row,]$Reads<=25000){
		metrics[row,]$Quality<- "shallow"
	}
}
#print NA values/feature
for (col in colnames(metrics)){
	message(col," has ",sum(is.na(metrics[,col]))," NAs")
}

#write to tables
write.table(metrics,"Output/Metrics/completeMetrics.txt",sep="\t",row.names = F,col.names = T,quote=F)
write.table(metrics$Library,"Output/TraningLibraries/trainingLibraries.txt",sep="\t",row.names = F,col.names = F,quote=F)
#remove NA rows and library feature
#metrics=na.exclude(metrics)
metrics=select(metrics,-c(Library))


for (row in 1:nrow(metrics)){
	if (metrics[row,]$Quality=="acceptable" | metrics[row,]$Quality=="a" ){
		metrics[row,]$Quality<- "good"
	} else if (metrics[row,]$Quality=="brdu" | metrics[row,]$Quality=="b" ){
		metrics[row,]$Quality="poor"
	} else if (metrics[row,]$Quality=="p"){
		metrics[row,]$Quality="poor"
	}
}
metrics$Quality <- droplevels.factor(metrics$Quality)
metrics %>% group_by(Quality)%>% summarize(n())


################## Data partitioning ##################
# create a list of 80% of the rows in the original dataset we can use for training


#nrow(na.omit(metrics))
metrics <- na.omit(metrics)
validation_index <- createDataPartition(metrics$Quality, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

#dataset=na.exclude(dataset)


################## Pre processing ###################

#Note imputation made performance on validation set worse even though there are no NAs
#Meaning norrmalizatiton was having a negative effect

#preProcess_missingdata_model <- preProcess(trainData, method='knnImpute')
#trainData <- predict(preProcess_missingdata_model, newdata = trainData)  #  !anyNA(trainData)
#preProcess_missingdata_model <- preProcess(testData, method='knnImpute')
#testData <- predict(preProcess_missingdata_model, newdata = testData)

#Note if we had any categorical variables we would need to one hot encode them
#But in this case we don't so our dataset is complete

################### Preliminary plotting ##################
#just  to see if any variables are obviously important

x <- trainData[,c(1:4,6:14)] #metrics
y <- trainData[,5] #classification
scales <- list(x=list(relation="free"), y=list(relation="free"))
#featurePlot(x=x, y=y, plot="box",scales =scales)
#featurePlot(x=x, y=y, plot="density",scales =scales,strip=strip.custom(par.strip.text=list(cex=.7)))
#featurePlot(x = metrics[, c("Background","Reads","complexity","spikiness","entropy")],y = metrics$Quality,plot = "pairs",auto.key = list())
#pairs(metrics[, c("complexity","spikiness","entropy","Background","Reads")],col=metrics$Quality)



################## Training models ##################

control <- trainControl(method="cv", number=10,classProbs = T)
metric <- "Sens"

#Linear (lda), non-linear(CART,kNN), advanced
models=c("lda","rpart","knn","svmRadial","rf","gbm")
#for (model in models){
#	set.seed(7)
#	fit <- train(Quality~., data=trainData, method=model, metric=metric, trControl=control)
#	assign(paste0("fit.",model),fit)
	#results=resamples(list(get(paste0("fit.",model))))
#}
fit.gbm = fit <- train(Quality~., data=trainData, method="gbm", metric=metric, trControl=control)

#results <- resamples(list(lda=fit.lda,rpart=fit.rpart,knn=fit.knn,svmRadial=fit.svmRadial,rf=fit.rf,gbm=fit.gbm))
#summary(results)
#dotplot(results)
#default_knn_mod = train(Quality ~ .,data = trainData,method = "knn",trControl = trainControl(method = "cv", number = 5),preProcess = c("center", "scale"),tuneGrid = expand.grid(k = seq(1, 101, by = 2)))
#fit.gbm.imp <- varImp(fit.gbm)
#plot(fit.rf.imp, main="Variable Importance with Random Forest")
plot(fit.gbm, main="Model Accuracies with Random Forest")

#Running best model on testData/test set
predictions <- predict(fit.gbm, na.omit(testData))
confusionMatrix(predictions, testData$Quality)



################## Feature selection ###################

set.seed(100)
options(warn=-1)
subsets <- c(1:13)
ctrl <- rfeControl(functions = rfFuncs,
				   method = "repeatedcv",
				   repeats = 5,
				   verbose = FALSE)

lmProfile <- rfe(x=trainData[, -5], y=trainData$Quality,
				 sizes = subsets,
				 rfeControl = ctrl)
lmProfile # the more variables the better apparerently

#Importance of features using Boruta package
fit.rf.Imp <- varImp(fit.rf, scale = FALSE)
plot(fit.rf.Imp)
boruta_output <- Boruta(Quality ~ ., data=na.omit(trainData), doTrace=0)
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
head(imps2[order(-imps2$meanImp), ])  # descending sort
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")


training_preprocess <-preProcess(select(trainData, - Quality),
								 method = c("center", "scale", "nzv", "pca"))
training_preprocess$rotation




