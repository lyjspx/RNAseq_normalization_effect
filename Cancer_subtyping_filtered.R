require('TCGAbiolinks')
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
## gene selection  ----------------------------------------------------
require(caret)
meanLargerThanAverage = apply(rawMRNA,1,sum) >  median(apply(rawMRNA,1,sum))
varLargerThanAverage = apply(rawMRNA,1,var) >  median(apply(rawMRNA,1,var))

rawCountMeetBoth <- rawMRNA[meanLargerThanAverage & varLargerThanAverage,]
dim(rawCountMeetBoth)

rawCountMeetBothGeneID <- 
  rawCountData@rowRanges$ensembl_gene_id[meanLargerThanAverage & varLargerThanAverage]

rawCountMeetBothCor = cor(t(rawCountMeetBoth))
rawCountMeetBothCor0.5 <- findCorrelation(rawCountMeetBothCor,cutoff = 0.5)
rm(rawCountMeetBothCor)
rawCountMeetBothCor0.5GeneID <- rawCountMeetBothGeneID[-rawCountMeetBothCor0.5]


##
meanLargerThanAverage = apply(fpkmMRNA,1,sum) > median(apply(fpkmMRNA,1,sum))
varLargerThanAverage = apply(fpkmMRNA,1,var) > median(apply(fpkmMRNA,1,var))

fpkmMeetBoth <- fpkmMRNA[meanLargerThanAverage & varLargerThanAverage,]
dim(fpkmMeetBoth)

fpkmMeetBothGeneID <- data@rowRanges$ensembl_gene_id[meanLargerThanAverage & varLargerThanAverage]

fpkmMeetBothCor = cor(t(fpkmMeetBoth))
fpkmMeetBothCor0.5 <- findCorrelation(fpkmMeetBothCor,cutoff = 0.5)
rm(fpkmMeetBothCor)
fpkmMeetBothGeneID[-fpkmMeetBothCor0.5]
fpkmMeetBothCor0.5GeneID <- fpkmMeetBothGeneID[-fpkmMeetBothCor0.5]


## edgeR filtering----------------------------------------------------------------
require('edgeR')
dgList <- DGEList(counts=rawMRNA,
                  genes=rawCountData@rowRanges$ensembl_gene_id)

countsPerMillion <- cpm(dgList)
countsPerMillion[1:5,1:5]

keep <- filterByExpr(countsPerMillion,min.count = 10,min.total.count = 0.7)

retainedEdgeR <- countsPerMillion[keep,]

dgList$genes[keep,]

dim(retainedEdgeR)
retainedEdgeR[1:5,1:5]
## pheno and geno matching -----------------------------------------------------------------
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
require('caret')
require('MLeval')
phenoFilteredOrdered$PAM50

#Sampling 10000 raw counts
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(rawMRNA)[,sample(1:56000,10000)])

#fpkmMeetBoth
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(fpkmMeetBoth[-fpkmMeetBothCor0.5,]))

#edgeR filtering
trainingSet <- cbind.data.frame(status=factor(phenoFilteredOrdered$PAM50),
                                t(retainedEdgeR))

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


trainingSet[1:5,1:5]
dim(trainingSet)

# Run algorithms using 5-fold cross validation
control <- trainControl(method="cv", number=5,savePredictions = TRUE, 
                        classProbs = TRUE, 
                        verboseIter = TRUE,allowParallel = TRUE)
metric <- c("Accuracy")

library(parallel)
library(doParallel)
require('mltest')
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

para_to_String <- function(fit.result,metrics){
  selected.model.index <- which(fit.result[metrics] == max(fit.result[metrics]))
  num.para <- which(colnames(fit.result) == metrics) - 1
  return(paste0(colnames(fit.result[1:num.para]),
                fit.result[selected.model.index,1:num.para],
                collapse = "_"))
}


set.seed(7) #fix random process to 'fixed'

# adaboost
fit.adaboost <- train(status~., data=trainingSet, method="adaboost", 
                      metric=metric, trControl=control,tuneLength=3)
fit.adaboost
# SVM
fit.svm <- train(status~., data=trainingSet, method="svmRadial", 
                 metric=metric, trControl=control,tuneLength=3)
fit.svm$results
mlTestSVM <- ml_test(fit.svm$pred$pred,fit.svm$pred$obs, output.as.table = TRUE)
mlTestSVM["Classifier"] <- "svmRadial"
mlTestSVM["Para"] <- para_to_String(fit.svm$results,"Accuracy")


fit.svmPoly <- train(status~., data=trainingSet, method="svmPoly", 
                     metric=metric, trControl=control,tuneLength=3)
fit.svmPoly$results
mlTestSVMPoly <- ml_test(fit.svmPoly$pred$pred,fit.svmPoly$pred$obs, output.as.table = TRUE)
mlTestSVMPoly["Classifier"] <- "svmPoly"
mlTestSVMPoly["Para"] <- para_to_String(fit.svmPoly$results,"Accuracy")

fit.rf <- train(status~., data=trainingSet, method="rf",
                metric=metric, trControl=control,tuneLength=3)
fit.rf$results
mlTestRF <-  ml_test(fit.rf$pred$pred,fit.rf$pred$obs, output.as.table = TRUE)
mlTestRF["Classifier"] <- "RF"
mlTestRF["Para"] <- para_to_String(fit.rf$results,"Accuracy")

fit.nerualnet <- train(status~., data=trainingSet, method="mlp",
                       metric=metric, trControl=control,tuneLength=3)
fit.nerualnet$finalModel
mlTestNN <- ml_test(fit.nerualnet$pred$pred,fit.nerualnet$pred$obs, output.as.table = TRUE)
mlTestNN["Classifier"] <- "NN"
mlTestNN["Para"] <- para_to_String(fit.nerualnet$results,"Accuracy")

fit.xgbTree <- train(status~., data=trainingSet, method="xgbTree",
                     metric=metric, trControl=control,tuneLength=1)
fit.xgbTree$results
mlTestXGboost <- ml_test(fit.xgbTree$pred$pred,fit.xgbTree$pred$obs, output.as.table = TRUE)
mlTestXGboost["Classifier"] <- "xgboost"
mlTestXGboost["Para"] <- para_to_String(fit.xgbTree$results,"Accuracy")

fit.glm <- train(status~., data=trainingSet, method="glm",
                 metric=metric, trControl=control,tuneLength=1)
fit.glm$results

fit.pls <- train(status~., data=trainingSet, method="pls",
                 metric=metric, trControl=control,tuneLength=30)
mlTestPLS <- ml_test(fit.pls$pred$pred,fit.pls$pred$obs, output.as.table = TRUE)
mlTestPLS["Classifier"] <- "PLS"

result <- rbind.data.frame(mlTestSVM,mlTestSVMPoly,
                           mlTestRF,mlTestXGboost,
                           mlTestNN,mlTestPLS)

stopCluster(cluster)

## TMM -------------------------------------------------
calcNormFactors(rawMRNA)

## Relative Log Expression -------------------------
calcFactorRLE(rawMRNA,p=0.9)

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

library(parallel)
library(doParallel)
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
set.seed(7) #fix random process to 'fixed'

# adaboost
fit.adaboost <- train(status~., data=trainingSet, method="adaboost", 
                      metric=metric, trControl=control,tuneLength=3)
# SVM
fit.svm <- train(status~., data=trainingSet, method="svmRadial", 
                 metric=metric, trControl=control,tuneLength=10)
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

library(parallel)
library(doParallel)
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

##  metastatis -------------------------------------------------------------
#fpkm intersection between fpkm0.5 and edgeR filtering
sum(fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,])
trainingSet <- cbind.data.frame(status=patientInfoMatched$classification_of_tumor,
                                t(fpkmMeetBoth[-fpkmMeetBothCor0.5,
                                ][fpkmMeetBothCor0.5GeneID %in% dgList$genes[keep,],]))
trainingSet <-  trainingSet[!is.na(trainingSet$status),]
dim(trainingSet)

# Run algorithms using 5-fold cross validation
control <- trainControl(method="cv", number=5,savePredictions = TRUE, 
                        classProbs = TRUE, 
                        verboseIter = TRUE,allowParallel = TRUE)
metric <- c("RMSE")

library(parallel)
library(doParallel)
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

fit.svmPoly <- train(status~., data=trainingSet, method="svmPoly", 
                 metric=metric, trControl=control,tuneLength=2)
fit.svmPoly$results


fit.rf <- train(status~., data=trainingSet, method="rf",
                metric=metric, trControl=control,tuneLength=6)
fit.rf$results

fit.nerualnet <- train(status~., data=trainingSet, method="mlp",
                       metric=metric, trControl=control,tuneLength=4)

fit.nerualnet$results


fit.xgbTree <- train(status~., data=trainingSet, method="xgbTree",
                     metric=metric, trControl=control,tuneLength=1)
fit.xgbTree$results


fit.glm <- train(status~., data=trainingSet, method="glm",
                     metric=metric, trControl=control,tuneLength=2)
fit.glm$results

fit.leapSeq <- train(status~., data=trainingSet, method="leapSeq",
                     metric=metric, trControl=control,tuneLength=2)
fit.leapSeq

fit.pls <- train(status~., data=trainingSet, method="pls",
                     metric=metric, trControl=control,tuneLength=1)
fit.pls

fit.ridge <- train(status~., data=trainingSet, method="ridge",
                   metric=metric, trControl=control,tuneLength=2)
fit.ridge

fit.lasso <- train(status~., data=trainingSet, method="lasso",
                   metric=metric, trControl=control,tuneLength=2)
fit.lasso

#warning: pretty slow
fit.blasso <- train(status~., data=trainingSet, method="blasso",
                   metric=metric, trControl=control,tuneLength=2)
fit.blasso

fit.bridge <- train(status~., data=trainingSet, method="bridge",
                    metric=metric, trControl=control,tuneLength=2)
fit.bridge

stopCluster(cluster)