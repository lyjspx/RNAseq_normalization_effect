if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("DESeq2")

BiocManager::install("XVector")

BiocManager::install("XML")

BiocManager::install("TCGAbiolinks")
BiocManager::install("edgeR")


require('TCGAbiolinks')

## FPKM ----------------------------------
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

GDCdownload(query = query,directory = '/home/rstudio/rdata/GDCdata')

data <- GDCprepare(query = query,directory = '/home/rstudio/rdata/GDCdata')


## raw Count --------------------------------------------
rawCountQuery <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - Counts")

GDCdownload(query = rawCountQuery,directory = '/home/rstudio/rdata/GDCdata')

rawCountData <- GDCprepare(query = rawCountQuery,directory = '/home/rstudio/rdata/GDCdata')

## FPKM-UQ ---------------------------------------------------

fpkmUqQuery <- GDCquery(project = "TCGA-BRCA",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - FPKM-UQ")

GDCdownload(query = fpkmUqQuery,directory = '/home/rstudio/rdata/GDCdata')

fpkmUqData <- GDCprepare(query = fpkmUqQuery,directory = '/home/rstudio/rdata/GDCdata')

## gene selection  ----------------------------------------------------
require(caret)
meanLargerThanAverage = apply(rawCountData@assays@data@listData$`HTSeq - Counts`,1,sum) >
  median(apply(rawCountData@assays@data@listData$`HTSeq - Counts`,1,sum))

varLargerThanAverage = apply(rawCountData@assays@data@listData$`HTSeq - Counts`,1,var) >
  median(apply(rawCountData@assays@data@listData$`HTSeq - Counts`,1,var))

rawCountMeetBoth <- 
  rawCountData@assays@data@listData$`HTSeq - Counts`[meanLargerThanAverage & varLargerThanAverage,]
dim(rawCountMeetBoth)

rawCountMeetBothGeneID <- rawCountData@rowRanges$ensembl_gene_id[meanLargerThanAverage & varLargerThanAverage]


rawCountMeetBothCor = cor(t(rawCountMeetBoth))
rawCountMeetBothCor0.5 <- findCorrelation(rawCountMeetBothCor,cutoff = 0.5)
rawCountMeetBothCor0.5GeneID <- rawCountMeetBothGeneID[-rawCountMeetBothCor0.5]
##
meanLargerThanAverage = apply(data@assays@data@listData$`HTSeq - FPKM`,1,sum) >
  median(apply(data@assays@data@listData$`HTSeq - FPKM`,1,sum))

varLargerThanAverage = apply(data@assays@data@listData$`HTSeq - FPKM`,1,var) >
  median(apply(data@assays@data@listData$`HTSeq - FPKM`,1,var))

fpkmMeetBoth <- 
  data@assays@data@listData$`HTSeq - FPKM`[meanLargerThanAverage & varLargerThanAverage,]
dim(fpkmMeetBoth)

fpkmMeetBothGeneID <- data@rowRanges$ensembl_gene_id[meanLargerThanAverage & varLargerThanAverage]

fpkmMeetBothCor = cor(t(fpkmMeetBoth))

fpkmMeetBothCor0.5 <- findCorrelation(fpkmMeetBothCor,cutoff = 0.5)
fpkmMeetBothGeneID[-fpkmMeetBothCor0.5]
fpkmMeetBothCor0.5GeneID <- fpkmMeetBothGeneID[-fpkmMeetBothCor0.5]
## edgeR filtering----------------------------------------------------------------
require('edgeR')
dgList <- DGEList(counts=rawCountData@assays@data@listData$`HTSeq - Counts`,
                  genes=rawCountData@rowRanges$ensembl_gene_id)

countsPerMillion <- cpm(dgList)
countsPerMillion[1:5,1:5]

keep <- filterByExpr(countsPerMillion,min.count = 10,min.total.count = 0.7)

retainedEdgeR <- countsPerMillion[keep,]

dgList$genes[keep,]

dim(retainedEdgeR)
retainedEdgeR[1:5,1:5]

## venn diagram ---------------------------------------------------------------
require("VennDiagram")
unionSet <- union(fpkmMeetBothCor0.5GeneID,rawCountMeetBothCor0.5GeneID)
intersectSet <- intersect(fpkmMeetBothCor0.5GeneID,rawCountMeetBothCor0.5GeneID)

intersectSet <- intersect(dgList$genes[keep,],rawCountMeetBothCor0.5GeneID)

venn.plot <- draw.triple.venn(area1 =  length(rawCountMeetBothCor0.5GeneID),
                              area2 =  length(fpkmMeetBothCor0.5GeneID),
                              area3 =  length(dgList$genes[keep,]),
                              n12 = length(intersect(fpkmMeetBothCor0.5GeneID,rawCountMeetBothCor0.5GeneID)),
                              n23 = length(intersect(fpkmMeetBothCor0.5GeneID,dgList$genes[keep,])),
                              n13 = length(intersect(rawCountMeetBothCor0.5GeneID,dgList$genes[keep,])),
                              n123 = length(intersect(intersect(fpkmMeetBothCor0.5GeneID,rawCountMeetBothCor0.5GeneID),
                                                      dgList$genes[keep,])),
                              c("Raw count", "FPKM", "edgeR filtering"))
grid.draw(venn.plot)

## prediction -------------------------------------------------------------
pheno <- read.csv('/home/rstudio/rdata/mmc2-1.csv',header = TRUE,stringsAsFactors = FALSE)
pheno[1:5,1:5]
colnames(pheno)


dim(query$results[[1]])
query$results[[1]][1:5,1:15]
query$results[[1]]$submitter_id

length(query$results[[1]]$cases)
length(query$results[[1]]$submitter_id)

length(unique(query$results[[1]]$cases))

length(unique(query$results[[1]]$cases.submitter_id))

length(intersect(pheno$mRNA ,query$results[[1]]$cases))

length(unique(query$results[[1]]$cases))

phenoAvaCase <- query$results[[1]]$cases %in% pheno$mRNA

