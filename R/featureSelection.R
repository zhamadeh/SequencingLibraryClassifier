
################## Feature selection ###################

#set.seed(100)
#options(warn=-1)
#subsets <- c(1:13)
#ctrl <- rfeControl(functions = rfFuncs,
#				   method = "repeatedcv",
#				   repeats = 5,
#				   verbose = FALSE)

#lmProfile <- rfe(x=trainData[, -5], y=trainData$Quality,
#				 sizes = subsets,
#				 rfeControl = ctrl)
#lmProfile # the more variables the better apparerently

#Importance of features using Boruta package
#fit.rf.Imp <- varImp(fit.rf, scale = FALSE)
#plot(fit.rf.Imp)
#boruta_output <- Boruta(Quality ~ ., data=na.omit(trainData), doTrace=0)
#boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
#roughFixMod <- TentativeRoughFix(boruta_output)
#boruta_signif <- getSelectedAttributes(roughFixMod)
#imps <- attStats(roughFixMod)
#imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
#head(imps2[order(-imps2$meanImp), ])  # descending sort
#plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")


#training_preprocess <-preProcess(select(trainData, - Quality),
#								 method = c("center", "scale", "nzv", "pca"))
#training_preprocess$rotation
