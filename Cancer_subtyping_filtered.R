library('TCGAbiolinks')
library('caret')
library('edgeR')
library('parallel')
library('doParallel')
library('MLeval')
library('mltest')

source('utils.R')
## data filtering --------------------------
pheno <- read.csv('/home/rstudio/rdata/mmc2-1.csv',header = TRUE,stringsAsFactors = FALSE)
pheno[1:5,1:5]
colnames(pheno)

dim(query$results[[1]])
query$results[[1]][1:5,1:15]
query$results[[1]]$submitter_id

length(query$results[[1]]$cases)
length(query$results[[1]]$submitter_id)

length(intersect(pheno$mRNA ,query$results[[1]]$cases))

length(unique(query$results[[1]]$cases))

phenoAvaCase <- query$results[[1]]$cases %in% pheno$mRNA

genoAvaCase <- pheno$mRNA %in% query$results[[1]]$cases

phenoFiltered <- pheno[genoAvaCase,]

patientInfo <- read.delim('/home/rstudio/rdata/clinical.tsv',sep = '\t',
                          header = TRUE,stringsAsFactors = FALSE)

patientInfo <- patientInfo[!duplicated(patientInfo$case_submitter_id),]
dim(patientInfo)
## FPKM ----------------------------------
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

#GDCdownload(query = query,directory = '/home/rstudio/rdata/GDCdata')
data <- GDCprepare(query = query,directory = '/home/rstudio/rdata/GDCdata')

data@assays@data@listData$`HTSeq - FPKM`[1:5,1:5]

fpkmMRNA <- data@assays@data@listData$`HTSeq - FPKM`[,phenoAvaCase]
## raw Count --------------------------------------------
rawCountQuery <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - Counts")

#GDCdownload(query = rawCountQuery,directory = '/home/rstudio/rdata/GDCdata')
rawCountData <- GDCprepare(query = rawCountQuery,directory = '/home/rstudio/rdata/GDCdata')

rawMRNA <- rawCountData@assays@data@listData$`HTSeq - Counts`[,phenoAvaCase]
## FPKM-UQ ---------------------------------------------------

fpkmUqQuery <- GDCquery(project = "TCGA-BRCA",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - FPKM-UQ")

#GDCdownload(query = fpkmUqQuery,directory = '/home/rstudio/rdata/GDCdata')

fpkmUqData <- GDCprepare(query = fpkmUqQuery,directory = '/home/rstudio/rdata/GDCdata')

fpkmUqMRNA <- fpkmUqData@assays@data@listData$`HTSeq - FPKM-UQ`[,phenoAvaCase]
## raw count gene selection  ---------------------------------------------------
raw.count.selected <- select.by.correlation(rawMRNA,
                                             rawCountData@rowRanges$ensembl_gene_id)
## fpkm gene selection  ---------------------------------------------------
fpkm.count.selected <- select.by.correlation(fpkmMRNA,
                                             rawCountData@rowRanges$ensembl_gene_id)
## fpkm-UQ gene selection  ---------------------------------------------------
fpkmUQ.count.selected <- select.by.correlation(fpkmUqMRNA,
                                             rawCountData@rowRanges$ensembl_gene_id)

## edgeR cpm filtering----------------------------------------------------------------

dgList <- DGEList(counts=rawMRNA,
                  genes=rawCountData@rowRanges$ensembl_gene_id)

countsPerMillion <- cpm(dgList)
countsPerMillion[1:5,1:5]
keep <- filterByExpr(countsPerMillion,min.count = 10,min.total.count = 0.7)
retainedEdgeR <- countsPerMillion[keep,]
dgList$genes[keep,]
dim(retainedEdgeR)
retainedEdgeR[1:5,1:5]

## TMM -------------------------------------------------
TMM.count <- DGEList(counts = rawMRNA)
TMM.count <- calcNormFactors(TMM.count,method = "TMM")
TMM.count <- cpm(TMM.count,log = F)
TMM.count[1:5,1:5]
dim(TMM.count)

## TMM selection -------------------------------------------------
TMM.count.selected <- select.by.correlation(TMM.count,
                                            rawCountData@rowRanges$ensembl_gene_id)

## Relative Log Expression RLE -----------------------------------
RLE.count <- DGEList(rawMRNA)
RLE.count <- calcNormFactors(RLE.count,method = "RLE")
RLE.count$counts[1:5,1:5]
RLE.count <- cpm(RLE.count,log = F)
RLE.count[1:5,1:5]
dim(RLE.count)

## RLE gene selection -----------------------------------
RLE.count.selected <- select.by.correlation(RLE.count,
                                            rawCountData@rowRanges$ensembl_gene_id)

## pheno and geno matching ---------------------------------------
query$results[[1]]$cases[phenoAvaCase]

phenoFilteredOrdered <- 
  phenoFiltered[match(query$results[[1]]$cases[phenoAvaCase],phenoFiltered$mRNA),]

patientInfoMatched <-  patientInfo[match(query$results[[1]]$cases.submitter_id[phenoAvaCase]
      ,patientInfo$case_submitter_id),]

#pheno and geno now matching
phenoFilteredOrdered$mRNA[1:5] == query$results[[1]]$cases[phenoAvaCase][1:5] 

phenoFilteredOrdered$Case.ID[1:10] == patientInfoMatched$case_submitter_id[1:10]

temp <- patientInfo[patientInfo$case_submitter_id %in% 
              query$results[[1]]$cases.submitter_id[phenoAvaCase],]


## Predicting with caret subtyping ------------------------------------------------------

phenoFilteredOrdered$PAM50

#raw counts
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(rawMRNA[sample(1:56602,10000,replace = FALSE),]))
#fpkm
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkmMRNA[sample(1:56602,10000,replace = FALSE),]))
#fpkm-Uq
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkmUqMRNA[sample(1:56602,10000,replace = FALSE),]))

#edgeR filtering
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(retainedEdgeR))
#raw counts selected
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(raw.count.selected$expression))

#fpkm selected
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkm.count.selected$expression))

#fpkm-UQ selected
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkmUQ.count.selected$expression))

#TMM
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(TMM.count[sample(1:56602,10000,replace = FALSE),]))
#TMM selected
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(TMM.count.selected$expression))
#RLE
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(RLE.count[sample(1:56602,10000,replace = FALSE),]))
#RLE selected
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(RLE.count.selected$expression))
trainingSet[1:5,1:5]
dim(trainingSet)

# Run algorithms using 5-fold cross validation
control <- trainControl(method="cv", number=5,savePredictions = TRUE, 
                        classProbs = TRUE, 
                        verboseIter = TRUE,allowParallel = TRUE)
metric <- c("Accuracy")
cluster <- makeCluster(detectCores() - 16) 
registerDoParallel(cluster)

set.seed(7) #fix random process to 'fixed'

# adaboost - unable apply in subtyping
# fit.adaboost <- train(status~., data=trainingSet, method="adaboost", 
#                       metric=metric, trControl=control,tuneLength=3)
# fit.adaboost
# SVM
fit.svm <- train(status~., data=trainingSet, method="svmRadial", 
                 metric=metric, trControl=control,tuneLength=3)
mlTestSVM <- ml_test(fit.svm$pred$pred,fit.svm$pred$obs, output.as.table = TRUE)
mlTestSVM["Classifier"] <- "svmRadial"
mlTestSVM["Para"] <- para_to_String(fit.svm$results,"Accuracy")


fit.svmPoly <- train(status~., data=trainingSet, method="svmPoly", 
                     metric=metric, trControl=control,tuneLength=3)
mlTestSVMPoly <- ml_test(fit.svmPoly$pred$pred,fit.svmPoly$pred$obs, output.as.table = TRUE)
mlTestSVMPoly["Classifier"] <- "svmPoly"
mlTestSVMPoly["Para"] <- para_to_String(fit.svmPoly$results,"Accuracy")

fit.rf <- train(status~., data=trainingSet, method="rf",
                metric=metric, trControl=control,tuneLength=3)
mlTestRF <-  ml_test(fit.rf$pred$pred,fit.rf$pred$obs, output.as.table = TRUE)
mlTestRF["Classifier"] <- "RF"
mlTestRF["Para"] <- para_to_String(fit.rf$results,"Accuracy")

fit.nerualnet <- train(status~., data=trainingSet, method="mlp",
                       metric=metric, trControl=control,tuneLength=3)
mlTestNN <- ml_test(fit.nerualnet$pred$pred,fit.nerualnet$pred$obs, output.as.table = TRUE)
mlTestNN["Classifier"] <- "NN"
mlTestNN["Para"] <- para_to_String(fit.nerualnet$results,"Accuracy")

fit.xgbTree <- train(status~., data=trainingSet, method="xgbTree",
                     metric=metric, trControl=control,tuneLength=1)
mlTestXGboost <- ml_test(fit.xgbTree$pred$pred,fit.xgbTree$pred$obs, output.as.table = TRUE)
mlTestXGboost["Classifier"] <- "xgboost"
mlTestXGboost["Para"] <- para_to_String(fit.xgbTree$results,"Accuracy")

# fit.glm <- train(status~., data=trainingSet, method="glm",
#                  metric=metric, trControl=control,tuneLength=1)
# fit.glm$results

fit.pls <- train(status~., data=trainingSet, method="pls",
                 metric=metric, trControl=control,tuneLength=30)
mlTestPLS <- ml_test(fit.pls$pred$pred,fit.pls$pred$obs, output.as.table = TRUE)
mlTestPLS["Classifier"] <- "PLS"
mlTestPLS["Para"] <- para_to_String(fit.pls$results,"Accuracy")

result <- rbind.data.frame(mlTestSVM,mlTestSVMPoly,
                           mlTestRF,mlTestXGboost,
                           mlTestNN,mlTestPLS)
write.csv(result,
          file = "/home/rstudio/rdata/Normalization_effect/results/subtyping/rle.selected.csv")

stopCluster(cluster)



## TumorPurity--------------------------------------

#fpkm intersection between fpkm0.5 and edgeR filtering
sum(fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,])
trainingSet <- cbind.data.frame(status=phenoFilteredOrdered$TumorPurity,
                                t(fpkmMeetBoth[-fpkmMeetBothCor0.5,][fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,],]))
trainingSet <-  trainingSet[!is.na(trainingSet$status),]

# Run algorithms using 5-fold cross validation
control <- trainControl(method="cv", number=5,savePredictions = TRUE, 
                        classProbs = FALSE, 
                        verboseIter = TRUE,allowParallel = TRUE)
metric <- c("RMSE")

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
set.seed(7) #fix random process to 'fixed'

# adaboost
fit.adaboost <- train(status~., data=trainingSet, method="adaboost", 
                      metric=metric, trControl=control,tuneLength=3)
# SVM
fit.svm <- train(status~., data=trainingSet, method="svmRadial", 
                 metric=metric, trControl=control,tuneLength=3)
fit.svm$results


fit.rf <- train(status~., data=trainingSet, method="rf",
                metric=metric, trControl=control,tuneLength=6)
fit.rf$results

fit.nerualnet <- train(status~., data=trainingSet, method="mlp",
                       metric=metric, trControl=control,tuneLength=6)

fit.nerualnet$result


fit.xgbTree <- train(status~., data=trainingSet, method="xgbTree",
                     metric=metric, trControl=control,tuneLength=1)
fit.xgbTree$results


fit.bridge <- train(status~., data=trainingSet, method="bridge",
                    metric=metric, trControl=control,tuneLength=1)
fit.bridge

fit.ridge <- train(status~., data=trainingSet, method="ridge",
                   metric=metric, trControl=control,tuneLength=1)
fit.ridge

fit.bridge <- train(status~., data=trainingSet, method="bridge",
                    metric=metric, trControl=control,tuneLength=1)
fit.bridge

fit.lasso <- train(status~., data=trainingSet, method="lasso",
                   metric=metric, trControl=control,tuneLength=1)
fit.lasso

stopCluster(cluster)
## Proliferation score--------------------------------------

#fpkm intersection between fpkm0.5 and edgeR filtering
sum(fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,])
trainingSet <- cbind.data.frame(status=phenoFilteredOrdered$ProliferationScore,
                                t(fpkmMeetBoth[-fpkmMeetBothCor0.5,
                                               ][fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,],]))
trainingSet <-  trainingSet[!is.na(trainingSet$status),]


# Run algorithms using 5-fold cross validation
control <- trainControl(method="cv", number=5,savePredictions = TRUE, 
                        classProbs = FALSE, 
                        verboseIter = TRUE,allowParallel = TRUE)
metric <- c("RMSE")
cluster <- makeCluster(detectCores() - 10) # convention to leave 1 core for OS
registerDoParallel(cluster)
set.seed(7) #fix random process to 'fixed'

# adaboost
fit.adaboost <- train(status~., data=trainingSet, method="adaboost", 
                      metric=metric, trControl=control,tuneLength=3)
# SVM
fit.svm <- train(status~., data=trainingSet, method="svmRadial", 
                 metric=metric, trControl=control,tuneLength=2)
fit.svm$results


fit.rf <- train(status~., data=trainingSet, method="rf",
                metric=metric, trControl=control,tuneLength=6)
fit.rf$results
fit.rf$modelType
fit.rf$pred
fit.rf$pred[fit.rf$pred$mtry==32,1] 
fit.rf$bestTune
fit.rf$finalModel$predicted
fit.rf$finalModel$mtry

fit.nerualnet <- train(status~., data=trainingSet, method="mlp",
                       metric=metric, trControl=control,tuneLength=4)

fit.nerualnet$results
fit.nerualnet$finalModel$fitted.values

fit.xgbTree <- train(status~., data=trainingSet, method="xgbTree",
                     metric=metric, trControl=control,tuneLength=2)
fit.xgbTree$results

fit.glm <- train(status~., data=trainingSet, method="glm",
                 metric=metric, trControl=control,tuneLength=2)
fit.glm$results

fit.knn <- train(status~., data=trainingSet, method="knn",
                 metric=metric, trControl=control,tuneLength=2)
fit.knn$results

fit.leapSeq <- train(status~., data=trainingSet, method="leapSeq",
                     metric=metric, trControl=control,tuneLength=20)
fit.leapSeq

fit.pls <- train(status~., data=trainingSet, method="pls",
                 metric=metric, trControl=control,tuneLength=10)
fit.pls

fit.bridge <- train(status~., data=trainingSet, method="bridge",
                    metric=metric, trControl=control,tuneLength=1)
fit.bridge

fit.ridge <- train(status~., data=trainingSet, method="ridge",
                   metric=metric, trControl=control,tuneLength=1)
fit.ridge

fit.bridge <- train(status~., data=trainingSet, method="bridge",
                    metric=metric, trControl=control,tuneLength=1)
fit.bridge

fit.lasso <- train(status~., data=trainingSet, method="lasso",
                   metric=metric, trControl=control,tuneLength=1)
fit.lasso

stopCluster(cluster)

## no more used - storage only ---------------------------------------------
#fpkm intersection between fpkm0.5 and edgeR filtering
sum(fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,])
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkmMeetBoth[-fpkmMeetBothCor0.5,
                                ][fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,],]))

#cpm intersection between fpkm0.5 and edgeR filtering
sum(dgList$genes[keep,] %in% fpkmMeetBothCor0.5GeneID)
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(retainedEdgeR[dgList$genes[keep,
                                ] %in% fpkmMeetBothCor0.5GeneID,]))